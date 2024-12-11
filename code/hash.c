#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "hash.h"
#include "blast.h"

// 哈希函数 (kmer -> hash value)
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

// 初始化哈希表
void init_hash_table(HashList* hash_table) {
    for (int i = 0; i < HASH_TABLE_SIZE; i++) {
        hash_table[i].head = NULL;
        hash_table[i].tail = NULL;
    }
}

// 建立哈希表
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

// 释放哈希表
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