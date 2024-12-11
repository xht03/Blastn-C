#include<stdlib.h>
#include<string.h>
#include<stdio.h>
#include<time.h>
#include "read.h"
#include "hash.h"
#include "blast.h"

#define QUERY_PATH "../seq/query.fasta"        // 查询序列的路径
#define SAMPLE_PATH "../seq/chr6.fasta"        // 样本序列的路径
#define OUTPUT_PATH "../output"                // 输出文件夹的路径
#define RESULT_NUM 20                          // 输出比对结果的数量
#define OVERLAP_THRESHOLD 0.9                  // 重叠阈值

HashList hash_table[HASH_TABLE_SIZE];          // 存放样本序列的 hash 表

int main()
{

    // ----------------- 读取样本序列名字列表 -----------------

    printf("开始读取\033[1m样本序列列表\033[0m\r");

    char** sample_names = read_seq_names(SAMPLE_PATH);      // 读取样本序列的名字
    int sample_num = count_seq(SAMPLE_PATH);                // 统计样本序列的数量

    printf("读取\033[32m样本序列列表\033[0m完成 样本序列数量: %d\n", sample_num);
    printf("\n");

    // ----------------- 读取查询序列 -----------------

    printf("开始读取\033[33m查询序列\033[0m\r");

    char* query = read_a_sequence(QUERY_PATH, "query");      // 读取查询序列
    int query_length = strlen(query);                       // 查询序列的长度

    printf("读取\033[32m查询序列\033[0m完成\n");
    printf("\n");

    // ----------------- 依次比对每个样本序列 -----------------

    printf("总共将执行 %d 个比对任务\n", sample_num);
    printf("\n");

    for (int i = 0; i < sample_num; i++) {

        // 初始化最终比对结果 (保存最高的 10 个比对结果)
        AlignmentResult top[RESULT_NUM];     
        for(int j = 0; j < RESULT_NUM; j++) {
            top[j].max_score = 0;
            top[j].start_pos1 = 0;
            top[j].end_pos1 = 0;
            top[j].start_pos2 = 0;
            top[j].end_pos2 = 0;
        }

        // ----------------- 读取样本序列 -----------------

        printf("-----\033[32m开始第 %d 个比对任务 %s \033[0m-----\n", i + 1, sample_names[i]);
        printf("1.正在读取序列 %s\r", sample_names[i]);
        fflush(stdout);

        char* sample = read_a_sequence(SAMPLE_PATH, sample_names[i]);   // 读取样本序列
        int sample_length = strlen(sample);                             // 样本序列的长度

        printf("1.读取序列 %s 完成\n", sample_names[i]);
        fflush(stdout);

        // ----------------- 建立样本序列的 hash 表 -----------------

        printf("2.开始建立 hash 表\r");
        fflush(stdout);

        init_hash_table(hash_table);                                    // 初始化 hash 表
        build_hash_table(hash_table, sample, sample_names[i]);          // 建立样本序列的 hash 表

        printf("2.建立 hash 表完成\n");
        fflush(stdout);

        // ----------------- 从查询序列中提取 kmer 并比对 -----------------

        printf("3.开始比对样本序列 %s\r", sample_names[i]);
        fflush(stdout);

        for (int k = 0; k <= query_length - KMER_LENGTH; k += 1) {
            // 从查询序列中提取 kmer
            Kmer kmer;
            strncpy(kmer.word, query + k, KMER_LENGTH);
            kmer.word[KMER_LENGTH] = '\0';
            kmer.position = k;

            // 计算 kmer 的哈希值
            unsigned int index = hash(kmer.word);

            // 在样本序列的 hash 表中查找
            HashNode* node = hash_table[index].head;

            while (node) {
                // 比对
                AlignmentResult result = extend_align(query, sample, query_length, sample_length, kmer.position, node->position);

                // 保存得分最高的 RESULT_NUM 个比对结果
                if (result.max_score >= SCORE_THRESHOLD) {
                    int flag = 0;                   // 用于标记是否有重叠
                    int overlap[RESULT_NUM] = {0};  // 用于标记是否与前面的比对结果有重叠

                    for (int u = 0; u < RESULT_NUM; u++) {
                        overlap[u] = is_overlap(result, top[u], OVERLAP_THRESHOLD);
                        flag |= overlap[u];
                    }

                    // 如果没有重叠则保存到 top 中
                    if(!flag) {
                        for (int u = 0; u < RESULT_NUM; u++) {
                            if (result.max_score > top[u].max_score) {
                                for (int v = RESULT_NUM - 1; v > u; v--) {
                                    top[v] = top[v - 1];
                                }
                                top[u] = result;
                                break;
                            }
                        }
                    }
                    // 如果高度重叠，则替换在相似结果中替换
                    else {
                        for (int u = 0; u < RESULT_NUM; u++) {
                            if (overlap[u]) {
                                if (result.max_score > top[u].max_score) {
                                    top[u] = result;
                                    break;
                                }
                            }
                        }
                    }
                }

                node = node->next;  // 下一个节点

            }

            printf("3.开始比对样本序列 %s : %.2f%%\r", sample_names[i], (k + 1) * 100.0 / (query_length - KMER_LENGTH));
            fflush(stdout);
        }

        printf("3.比对完成\n");

        // ----------------- 输出比对结果 -----------------

        char output_filename[100];
        sprintf(output_filename, "%s/%s", OUTPUT_PATH, sample_names[i]);

        FILE* output_file = fopen(output_filename, "w");
        if (output_file) {
            for (int j = 0; j < RESULT_NUM; j++) {
                if (top[j].max_score > 0) {
                    fprintf(output_file, "比对结果 %d:\n", j + 1);
                    fprintf(output_file, "得分: %lld\n", top[j].max_score);
                    fprintf(output_file, "查询序列位置: [%d, %d]\n", top[j].start_pos1, top[j].end_pos1);
                    fprintf(output_file, "样本序列位置: [%d, %d]\n", top[j].start_pos2, top[j].end_pos2);
                    fprintf(output_file, "\n");
                }
            }
            fclose(output_file);
        } 
        else {
            perror("无法创建输出文件");
        }

        // ----------------- 释放内存 -----------------

        free(sample);                       // 释放样本序列
        release_hash_table(hash_table);     // 释放 hash 表

        printf("-----\033[32m第 %d 个比对任务 %s 完成\033[0m-----\n", i + 1, sample_names[i]);
        printf("\n");
    }

    printf("\033[32m所有任务完成\033[0m\n");

    return 0;
    
}