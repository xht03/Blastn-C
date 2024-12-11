#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

typedef struct SeqNode
{
    char* name;                 // 序列名
    char* sequence;             // 序列内容
    size_t length;              // 序列长度
    struct SeqNode* next;       // 下一个节点
} SeqNode;

char* check_name(char* line);                           // 从描述行中提取序列名
int count_seq(char* filename);                          // 统计文件中的序列数量
char** read_seq_names(char* filename);                  // 读取序列名字列表
char* read_a_sequence(char* filename, char* seqname);   // 读取一个序列
SeqNode* read_all_sequences(char* filename);            // 读取所有序列，并存储到链表中