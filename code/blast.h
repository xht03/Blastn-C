#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define KMER_LENGTH 11              // kmer 长度
#define KMER_STEP 100               // kmer 提取步长 (从样本序列中提取)
#define EXTEND_LENGTH 100           // 扩展比对的长度
#define EXTEND_THRESHOLD 350        // 扩展阈值 (当扩展切片得分低于此值时停止扩展)
#define SCORE_THRESHOLD 1000        // 得分阈值 (低于此值的比对结果将被忽略)

#define MATCH_SCORE 5           // 匹配得分
#define MISMATCH_SCORE -3       // 错配得分
#define GAP_PENALTY -2          // 空位惩罚

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MAX3(a, b, c) (MAX(MAX((a), (b)), (c)))
#define MAX4(a, b, c, d) (MAX(MAX3((a), (b), (c)), (d)))

typedef struct 
{
    char word[KMER_LENGTH + 1];  // kmer
    int position;                // 在序列上的位置
} Kmer;

typedef struct {
    int score;      // 得分
    int direction;  // 0: none, 1: diagonal, 2: up, 3: left
} Cell;

typedef struct {
    long long max_score;
    unsigned int start_pos1;
    unsigned int end_pos1;
    unsigned int start_pos2;
    unsigned int end_pos2;
} AlignmentResult;

AlignmentResult needleman_wunsch(char* seq1, char* seq2);
AlignmentResult extend_align(char* query, char* sample, int query_length, int sample_length, int query_start, int sample_start);
int is_overlap(AlignmentResult result1, AlignmentResult result2, float threshold);