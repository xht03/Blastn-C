#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define HASH_TABLE_SIZE 1048576     // 哈希表大小 4^11

typedef struct HashNode {
    int position;               // 序列上的位置
    char* sequence;             // 在那个序列
    struct HashNode* next;      // 下一个节点
} HashNode;

typedef struct {
    HashNode* head;         // 头节点
    HashNode* tail;         // 尾节点
} HashList;

unsigned int hash(char* kmer);
void init_hash_table(HashList* hash_table);
void build_hash_table(HashList* hash_table, char* sequence, char* seqname);
void release_hash_table(HashList* hash_table);



