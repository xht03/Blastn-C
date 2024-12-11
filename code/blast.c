#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "blast.h"


// 半全局比对 (与 Needleman-Wunsch 类似)
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
    for (int i = 1; i <= len1; i++) {
        for (int j = 1; j <= len2; j++) {
            // 计算得分
            int match = score_matrix[i-1][j-1].score + (seq1[i-1] == seq2[j-1] ? MATCH_SCORE : MISMATCH_SCORE);
            int delet = score_matrix[i-1][j].score + GAP_PENALTY;
            int insert = score_matrix[i][j-1].score + GAP_PENALTY;
            int score = MAX3(match, delet, insert);

            score_matrix[i][j].score = score;
        }
    }

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

// 扩展比对 (从 kmer 的位置开始左右扩展)
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
    
    //* 左侧扩展比对
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

    //* 右侧扩展比对
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

// 检测两个比对结果是否有重叠
int is_overlap(AlignmentResult result1, AlignmentResult result2, float threshold) {
    int overlap_start = MAX(result1.start_pos1, result2.start_pos1);
    int overlap_end = MIN(result1.end_pos1, result2.end_pos1);
    int overlap_length = overlap_end - overlap_start + 1;
    int result1_length = result1.end_pos1 - result1.start_pos1 + 1;
    int result2_length = result2.end_pos1 - result2.start_pos1 + 1;

    if (overlap_length > 0 && (overlap_length > threshold * result1_length && overlap_length > threshold * result2_length)) {
        return 1;
    }
    return 0;
}