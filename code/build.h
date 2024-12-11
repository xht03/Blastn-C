#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define KMER_LENGTH 11              // kmer 长度
#define KMER_STEP 100               // kmer 提取步长 (从样本序列中提取)
#define HASH_TABLE_SIZE 1048576     // 哈希表大小 4^11
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


typedef struct SeqNode
{
    char* name;           // 序列名
    char* sequence;       // 序列内容
    size_t length;        // 序列长度
    struct SeqNode* next;       // 下一个节点
} SeqNode;


typedef struct HashNode {
    int position;               // 序列上的位置
    char* sequence;             // 在那个序列
    struct HashNode* next;      // 下一个节点
} HashNode;

typedef struct {
    HashNode* head;         // 头节点
    HashNode* tail;         // 尾节点
} HashList;


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


// ----------------- 从FASTA文件中读取序列 -----------------


char* check_name(char* line) {
    char* name = NULL;
    char* space_pos = strchr(line, ' ');
    char* enter_pos = strchr(line, '\n');
    if (space_pos) {
        *space_pos = '\0';
    }
    if (enter_pos) {
        *enter_pos = '\0';
    }
    name = strdup(line + 1);  // 跳过 '>'
    return name;
}


int count_seq(char* filename) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        perror("无法打开文件");
        return -1;
    }

    int count = 0;
    char line[1024];

    while (fgets(line, sizeof(line), file)) {
        if (line[0] == '>') {
            count++;
        }
    }

    fclose(file);
    return count;
}


char** read_seq_names(char* filename) {
    // 打开文件
    FILE* file = fopen(filename, "r");
    if (!file) {
        perror("无法打开文件");
        return NULL;
    }

    char** name_list = NULL;    // 名字列表
    size_t list_size = 0;       // 名字数量
    
    char line[1024];            // 行缓冲区

    while (fgets(line, sizeof(line), file)) {

        // 找到新的序列 (描述行)
        if (line[0] == '>') {
            char* name = check_name(line);
            if (!name) {
                continue;
            }
            name_list = (char**)realloc(name_list, (list_size + 1) * sizeof(char*));
            name_list[list_size++] = name;
        }
    }

    fclose(file);

    return name_list;
}



char* read_a_sequence(char* filename, char* seqname) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        perror("无法打开文件");
        return NULL;
    }

    // 分配初始缓冲区
    size_t buffer_size = 1024;
    char* sequence = (char*)malloc(buffer_size);
    if (!sequence) {
        perror("内存分配失败");
        fclose(file);
        return NULL;
    }

    size_t length = 0;
    char line[1024];
    int found = 0;

    // 读取文件行
    while (fgets(line, sizeof(line), file)) {
        // 检查描述行
        if (line[0] == '>') {
            // 如果找到目标序列，停止读取
            if (found) {
                break;
            }
            // 检查描述行是否包含目标序列名
            if (strstr(line, seqname)) {
                found = 1;
            }
            continue;
        }

        // 如果还没有找到目标序列，继续读取
        if (!found) {
            continue;
        }

        // 去掉行末的换行符
        line[strcspn(line, "\n")] = '\0';

        // 检查缓冲区是否足够大 (不足则翻倍)
        size_t line_length = strlen(line);
        if (length + line_length >= buffer_size) {
            buffer_size *= 2;
            sequence = (char*)realloc(sequence, buffer_size);
            if (!sequence) {
                perror("内存重新分配失败");
                fclose(file);
                return NULL;
            }
        }

        // 追加行到序列
        strcpy(sequence + length, line);
        length += line_length;
    }

    fclose(file);

    // 如果没有找到目标序列，释放内存并返回NULL
    if (!found) {
        free(sequence);
        return NULL;
    }

    return sequence;
}


SeqNode* read_all_sequences(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        perror("无法打开文件");
        return NULL;
    }

    SeqNode* head = NULL;
    SeqNode* tail = NULL;
    char line[1024];
    size_t buffer_size = 1024;

    int cnt = 0;

    while (fgets(line, sizeof(line), file)) {

        cnt++;

        // 如果遇到新的描述行
        if (line[0] == '>') {
            // 创建新的序列节点
            SeqNode* node = (SeqNode*)malloc(sizeof(SeqNode));
            buffer_size = 1024;
            node->name = check_name(line);
            node->sequence = (char*)malloc(buffer_size);
            memset(node->sequence, 0, buffer_size);
            node->length = 0;
            node->next = NULL;

            if(head == NULL) {
                head = node;
                tail = node;
            } else {
                tail->next = node;
                tail = node;
            }
            
        }
        else {
            // 去掉行末的换行符
            line[strcspn(line, "\n")] = '\0';

            // 计算序列长度
            size_t line_length = strlen(line);
            if(tail->length + line_length >= buffer_size) {
                buffer_size *= 2;
                tail->sequence = (char*)realloc(tail->sequence, buffer_size);
            }
            // 追加行到序列
            strcpy(tail->sequence + tail->length, line);
            tail->length += line_length;
        }
    }

    fclose(file);
    return head;
}


// ----------------- 建立样本序列的 hash 表 -----------------


unsigned int hash(char* kmer) {
    unsigned int hash_value = 0;
    for (int i = 0; i < KMER_LENGTH; ++i) {
        hash_value <<= 2;
        switch (kmer[i]) {
            case 'A': hash_value |= 0; break;   // A -> 00
            case 'C': hash_value |= 1; break;   // C -> 01
            case 'G': hash_value |= 2; break;   // G -> 10
            case 'T': hash_value |= 3; break;   // T -> 11
        }
    }
    return hash_value % HASH_TABLE_SIZE;
}


void init_hash_table(HashList* hash_table) {
    for (int i = 0; i < HASH_TABLE_SIZE; i++) {
        hash_table[i].head = NULL;
        hash_table[i].tail = NULL;
    }
}


void build_hash_table(HashList* hash_table, char* sequence, char* seqname) {
    int length = strlen(sequence);
    for (int i = 0; i <= length - KMER_LENGTH; i += KMER_STEP) {
        // 从序列中提取 kmer
        char kmer[KMER_LENGTH + 1];
        strncpy(kmer, sequence + i, KMER_LENGTH);
        kmer[KMER_LENGTH] = '\0';

        // 创建新节点
        HashNode* node = (HashNode*)malloc(sizeof(HashNode));
        node->position = i;         // 从 0 开始
        node->sequence = seqname;
        node->next = NULL;

        // 计算 kmer 的哈希值
        unsigned int index = hash(kmer);

        // 链接到哈希表上
        if (hash_table[index].head == NULL) {
            hash_table[index].head = node;
            hash_table[index].tail = node;
        } else {
            hash_table[index].tail->next = node;
            hash_table[index].tail = node;
        }
    }
}


void release_hash_table(HashList* hash_table) {
    for (int i = 0; i < HASH_TABLE_SIZE; i++) {
        HashNode* node = hash_table[i].head;
        while (node) {
            HashNode* next = node->next;
            free(node);
            node = next;
        }
    }
}


// ----------------- 序列比对 -----------------


AlignmentResult needleman_wunsch(char* seq1, char* seq2) {
    int len1 = strlen(seq1);
    int len2 = strlen(seq2);

    // 创建得分矩阵
    Cell score_matrix[len1 + 1][len2 + 1];

    // 初始化得分矩阵
    for (int i = 0; i <= len1; i++) {
        score_matrix[i][0].score = i * GAP_PENALTY;
        score_matrix[i][0].direction = 2; // up
    }
    for (int j = 0; j <= len2; j++) {
        score_matrix[0][j].score = j * GAP_PENALTY;
        score_matrix[0][j].direction = 3; // left
    }

    // 初始化比对结果
    AlignmentResult result;
    result.max_score = 0;
    result.start_pos1 = 0;
    result.start_pos2 = 0;
    result.end_pos1 = 0;
    result.end_pos2 = 0;

    // 填充得分矩阵
    int max_i = 0, max_j = 0;
    for (int i = 1; i <= len1; i++) {
        for (int j = 1; j <= len2; j++) {
            // 计算得分
            int match = score_matrix[i-1][j-1].score + (seq1[i-1] == seq2[j-1] ? MATCH_SCORE : MISMATCH_SCORE);
            int delet = score_matrix[i-1][j].score + GAP_PENALTY;
            int insert = score_matrix[i][j-1].score + GAP_PENALTY;
            int score = MAX3(match, delet, insert);

            score_matrix[i][j].score = score;

            // 记录回溯方向
            // if (score == match) {
            //     score_matrix[i][j].direction = 1; // diagonal
            // } else if (score == delete) {
            //     score_matrix[i][j].direction = 2; // up
            // } else {
            //     score_matrix[i][j].direction = 3; // left
            // }
        }
    }

    // for debug
    // printf("得分矩阵如下：\n");
    // for (int i = 0; i <= len1; i++) {
    //     for (int j = 0; j <= len2; j++) {
    //         printf("%d ", score_matrix[i][j].score);
    //     }
    //     printf("\n");
    // }

    // 找到最大得分
    for (int i = 0; i <= len1; i++) {
        if(score_matrix[i][len2].score > result.max_score) {
            result.max_score = score_matrix[i][len2].score;
            result.end_pos1 = i;
            result.end_pos2 = len2;
        }
    }
    for (int j = 0; j <= len2; j++) {
        if(score_matrix[len1][j].score > result.max_score) {
            result.max_score = score_matrix[len1][j].score;
            result.end_pos1 = len1;
            result.end_pos2 = j;
        }
    }

    return result;
}


AlignmentResult extend_align(char* query, char* sample, int query_length, int sample_length, int query_start, int sample_start) {
    int query_end = query_start + KMER_LENGTH;      // 指向 query 中比对结束位置的下一个
    int sample_end = sample_start + KMER_LENGTH;    // 指向 sample 中比对结束位置的下一个
    
    // 初始化比对结果
    AlignmentResult result = {0, 0, 0, 0, 0};
    result.max_score = MATCH_SCORE * KMER_LENGTH;
    result.start_pos1 = query_start;
    result.start_pos2 = sample_start;
    result.end_pos1 = query_end;
    result.end_pos2 = sample_end;
    

    // 左侧扩展比对
    while(1) {
        // 检测两个左侧可切片的长度
        int left_query_length = (query_start > EXTEND_LENGTH ? EXTEND_LENGTH : query_start);
        int left_sample_length = (sample_start > EXTEND_LENGTH ? EXTEND_LENGTH : sample_start);
        if(left_query_length == 0 || left_sample_length == 0) {
            break;
        }

        // 获取 query 序列的左侧切片 (逆序)
        char* query_left = (char*)malloc(left_query_length + 1);
        for (int i = 0; i < left_query_length; i++) {
            query_left[i] = query[query_start - 1 - i];
        }
        query_left[left_query_length] = '\0';

        // 获取 sample 序列的左侧切片 (逆序)
        char* sample_left = (char*)malloc(left_sample_length + 1);
        for (int i = 0; i < left_sample_length; i++) {
            sample_left[i] = sample[sample_start - 1 - i];
        }
        sample_left[left_sample_length] = '\0'; 

        // 扩展比对
        AlignmentResult left_result = needleman_wunsch(query_left, sample_left);
        if (left_result.max_score < EXTEND_THRESHOLD || left_result.end_pos1 == 0 || left_result.end_pos2 == 0) {
            free(query_left);
            free(sample_left);
            break;
        }
        else {
            result.start_pos1 = query_start - left_result.end_pos1;
            result.start_pos2 = sample_start - left_result.end_pos2;
            result.max_score += left_result.max_score;

            query_start -= left_result.end_pos1;
            sample_start -= left_result.end_pos2;

            free(query_left);
            free(sample_left);
        }
    }

    // 右侧扩展比对
    while(1) {
        // 检测两个右侧可切片的长度
        int right_query_length = (query_length - query_end > EXTEND_LENGTH ? EXTEND_LENGTH : query_length - query_end);
        int right_sample_length = (sample_length - sample_end > EXTEND_LENGTH ? EXTEND_LENGTH : sample_length - sample_end);
        if(right_query_length == 0 || right_sample_length == 0) {
            break;
        }

        // 获取 query 序列的右侧切片
        char* query_right = (char*)malloc(right_query_length + 1);
        for (int i = 0; i < right_query_length; i++) {
            query_right[i] = query[query_end + i];
        }
        query_right[right_query_length] = '\0';

        // 获取 sample 序列的右侧切片
        char* sample_right = (char*)malloc(right_sample_length + 1);
        for (int i = 0; i < right_sample_length; i++) {
            sample_right[i] = sample[sample_end + i];
        }
        sample_right[right_sample_length] = '\0'; 

        // 扩展比对
        AlignmentResult right_result = needleman_wunsch(query_right, sample_right);
        if (right_result.max_score < EXTEND_THRESHOLD || right_result.end_pos1 == 0 || right_result.end_pos2 == 0) {
            free(query_right);
            free(sample_right);
            break;
        }
        else {
            result.end_pos1 = query_end + right_result.end_pos1;
            result.end_pos2 = sample_end + right_result.end_pos2;
            result.max_score += right_result.max_score;

            query_end += right_result.end_pos1;
            sample_end += right_result.end_pos2;

            free(query_right);
            free(sample_right);
        }
    }

    return result;
}








