package abra;

import htsjdk.samtools.*;
import htsjdk.samtools.util.Log;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicLong;

public class MultiThreadUnmappedReader {
    private static final Log log = Log.getInstance(MultiThreadUnmappedReader.class);
    private final File bamFile;
    private final File baiFile;
    private final int threadNum;
    private final AtomicLong totalUnmapped = new AtomicLong(0);
    private final AtomicLong totalProcessed = new AtomicLong(0);

    public MultiThreadUnmappedReader(File bamFile, File baiFile, int threadNum) {
        this.bamFile = bamFile;
        this.baiFile = baiFile;
        this.threadNum = threadNum;
    }

    /**
     * 多线程读取所有 unmapped reads
     * 策略：
     * 1. 按 contig 分块并行读取 mapped reads（快速跳过）
     * 2. 单独读取 unmapped reads（使用 queryUnmapped 专用接口）
     */
    public List<SAMRecord> readAllUnmapped() throws IOException, InterruptedException, ExecutionException {
        SamReader headerReader = null;
        try {
            // 检查索引文件
            if (baiFile == null || !baiFile.exists()) {
                log.warn("⚠️  未找到索引文件 " + baiFile + "，无法使用多线程，降级为单线程模式");
                return readAllUnmappedSingleThread();
            }
            
            // 创建 SamReaderFactory（只用于读取 header）
            SamReaderFactory factory = SamReaderFactory.makeDefault()
                    .validationStringency(ValidationStringency.SILENT);
            
            headerReader = factory.open(bamFile);
            
            // 获取参考序列信息
            SAMFileHeader header = headerReader.getFileHeader();
            List<SAMSequenceRecord> contigs = header.getSequenceDictionary().getSequences();
            
            if (contigs.isEmpty()) {
                log.warn("⚠️  BAM 文件中没有参考序列信息");
                return readAllUnmappedSingleThread();
            }
            
            log.info(String.format("📊 发现 %d 个参考序列，开始多线程处理...", contigs.size()));
            
            // 创建线程池
            ExecutorService executor = Executors.newFixedThreadPool(threadNum);
            List<Future<ContigResult>> futures = new ArrayList<>();
            
            // 提交任务：每个 contig 一个线程（每个线程打开自己的 reader）
            for (SAMSequenceRecord contig : contigs) {
                final String contigName = contig.getSequenceName();
                final int contigLength = contig.getSequenceLength();
                
                futures.add(executor.submit(() -> {
                    return processContig(contigName, contigLength);
                }));
            }
            
            // 提交 unmapped reads 读取任务（独立线程）
            Future<ContigResult> unmappedFuture = executor.submit(() -> {
                return processUnmappedReads();
            });
            futures.add(unmappedFuture);
            
            // 收集结果
            List<SAMRecord> allUnmapped = new CopyOnWriteArrayList<>();
            for (Future<ContigResult> future : futures) {
                ContigResult result = future.get();
                totalProcessed.addAndGet(result.processed);
                allUnmapped.addAll(result.unmapped);
            }
            
            executor.shutdown();
            executor.awaitTermination(1, TimeUnit.HOURS);
            
            totalUnmapped.set(allUnmapped.size());
            log.info(String.format("✅ 多线程完成！共处理 %d 条 reads，找到 %d 条 unmapped reads", 
                    totalProcessed.get(), allUnmapped.size()));
            
            return allUnmapped;
            
        } catch (Exception e) {
            log.error("多线程读取失败，降级为单线程模式", e);
            return readAllUnmappedSingleThread();
        } finally {
            if (headerReader != null) headerReader.close();
        }
    }
    
    /**
     * 处理单个染色体（每个线程创建独立的 reader）
     */
    private ContigResult processContig(String contigName, int contigLength) {
        ContigResult result = new ContigResult();
        result.contigName = contigName;
        
        SamReader reader = null;
        try {
            // 每个线程创建独立的 reader
            SamReaderFactory factory = SamReaderFactory.makeDefault()
                    .validationStringency(ValidationStringency.SILENT);
            reader = factory.open(SamInputResource.of(bamFile).index(baiFile));
            
            try (SAMRecordIterator iter = reader.query(contigName, 0, contigLength, false)) {
                while (iter.hasNext()) {
                    SAMRecord rec = iter.next();
                    result.processed++;
                    
                    // 检查是否为 unmapped read
                    if (rec.getReadUnmappedFlag()) {
                        result.unmapped.add(rec);
                    }
                }
            }
        } catch (Exception e) {
            log.error("处理染色体 " + contigName + " 失败", e);
        } finally {
            try {
                if (reader != null) reader.close();
            } catch (Exception e) {
                log.error("关闭 reader 失败", e);
            }
        }
        
        return result;
    }
    
    /**
     * 专门处理 unmapped reads（每个线程创建独立的 reader）
     */
    private ContigResult processUnmappedReads() {
        ContigResult result = new ContigResult();
        result.contigName = "unmapped";
        
        SamReader reader = null;
        try {
            // 每个线程创建独立的 reader
            SamReaderFactory factory = SamReaderFactory.makeDefault()
                    .validationStringency(ValidationStringency.SILENT);
            reader = factory.open(SamInputResource.of(bamFile).index(baiFile));
            
            try (SAMRecordIterator iter = reader.queryUnmapped()) {
                while (iter.hasNext()) {
                    SAMRecord rec = iter.next();
                    result.processed++;
                    result.unmapped.add(rec);
                }
            }
        } catch (Exception e) {
            log.error("读取 unmapped reads 失败", e);
        } finally {
            try {
                if (reader != null) reader.close();
            } catch (Exception e) {
                log.error("关闭 reader 失败", e);
            }
        }
        
        return result;
    }
    
    /**
     * 单线程降级方案
     */
    private List<SAMRecord> readAllUnmappedSingleThread() throws IOException {
        SamReader reader = null;
        SAMRecordIterator iter = null;
        try {
            SamReaderFactory factory = SamReaderFactory.makeDefault()
                    .validationStringency(ValidationStringency.SILENT);
            
            reader = factory.open(bamFile);
            
            log.info("开始单线程遍历 BAM 文件...");
            
            iter = reader.iterator();
            List<SAMRecord> result = new ArrayList<>();
            long count = 0;
            long total = 0;
            
            while (iter.hasNext()) {
                SAMRecord rec = iter.next();
                total++;
                
                if (rec.getReadUnmappedFlag()) {
                    result.add(rec);
                    count++;
                    
                    if (count % 10000 == 0) {
                        log.info(String.format("已处理 %d 条 reads，找到 %d 条 unmapped", total, count));
                    }
                }
            }
            
            log.info(String.format("✅ 单线程完成！共扫描 %d 条 reads，找到 %d 条 unmapped reads", total, result.size()));
            return result;
            
        } finally {
            try {
                if (iter != null) iter.close();
                if (reader != null) reader.close();
            } catch (IOException e) {
                log.error("关闭 reader 失败", e);
            }
        }
    }
    
    /**
     * 染色体处理结果
     */
    private static class ContigResult {
        String contigName;
        long processed = 0;
        List<SAMRecord> unmapped = new ArrayList<>();
    }

    /**
     * 写入 unmapped reads 到输出 BAM
     */
    public void writeUnmappedBam(List<SAMRecord> reads, File outFile) {
        if (reads.isEmpty()) {
            log.info("没有 unmapped reads 需要写入");
            return;
        }
        
        try (SamReader templateReader = SamReaderFactory.makeDefault().open(bamFile)) {
            SAMFileHeader header = templateReader.getFileHeader();
            
            try (SAMFileWriter writer = new SAMFileWriterFactory()
                    .setCreateIndex(true)
                    .setCreateMd5File(false)
                    .makeBAMWriter(header, true, outFile)) {
                
                for (SAMRecord rec : reads) {
                    writer.addAlignment(rec);
                }
            }
            log.info(String.format("✅ 已写入 %d 条 unmapped reads 到 %s", reads.size(), outFile));
            
        } catch (Exception e) {
            throw new RuntimeException("写入失败", e);
        }
    }

    // ============ 主方法 ============
    public static void main(String[] args) {
        if (args.length < 2) {
            System.err.println("用法: java MultiThreadUnmappedReader <input.bam> <output.bam> [线程数]");
            System.err.println("示例: java MultiThreadUnmappedReader input.bam output.bam 8");
            System.err.println("      java MultiThreadUnmappedReader input.bam output.bam (默认4线程)");
            System.exit(1);
        }

        File bamFile = new File(args[0]);
        File outFile = new File(args[1]);
        int threads = args.length >= 3 ? Integer.parseInt(args[2]) : 4;
        File baiFile = new File(args[0] + ".bai");

        if (!bamFile.exists()) {
            System.err.println("❌ BAM 文件不存在: " + bamFile);
            System.exit(1);
        }

        System.out.println("========================================");
        System.out.println("  多线程 BAM Unmapped Reads 提取工具");
        System.out.println("========================================");
        System.out.println("输入文件: " + bamFile);
        System.out.println("输出文件: " + outFile);
        System.out.println("线程数:   " + threads);
        System.out.println("索引文件: " + (baiFile.exists() ? "✅ 存在" : "❌ 不存在"));
        System.out.println("========================================");
        System.out.println();

        MultiThreadUnmappedReader reader = new MultiThreadUnmappedReader(bamFile, baiFile, threads);
        try {
            long start = System.currentTimeMillis();
            
            // 读取所有 unmapped reads
            List<SAMRecord> unmapped = reader.readAllUnmapped();
            
            // 写入结果
            reader.writeUnmappedBam(unmapped, outFile);
            
            long elapsed = (System.currentTimeMillis() - start) / 1000;
            System.out.println();
            System.out.println("========================================");
            System.out.printf("🎉 完成! 提取 %d 条 unmapped reads, 耗时 %d 秒%n", 
                            unmapped.size(), elapsed);
            System.out.println("========================================");
            
        } catch (Exception e) {
            log.error("执行失败", e);
            e.printStackTrace();
            System.exit(1);
        }
    }
}
