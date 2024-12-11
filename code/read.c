#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "read.h"


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


SeqNode* read_all_sequences(char* filename) {
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