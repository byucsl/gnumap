#ifndef ALIGN_H
#define ALIGN_H

enum ScoringType {
    BLOSUM62,
    DEFAULT
};

class Align
{
public:
    Align();

    double When();
    void initScoreMatrix(ScoringType type = DEFAULT);
    void loadMatchSequence(char* seq, int len);
    void loadRead(float* seq, int len);
    void loadSequences(char* seqs, int seqLen, int numSeq);
    void setGapPenalty(float penalty);
    bool setGPU(int num);
    void computeAlignments();
    void computeDenominator(float& denom);
    float* getScores();


private:
    char* seqA;
    int seqALen;
    char* seqB;
    float* seqBPtr;
    int seqBLen;
    int numSeq;
    float* scores;
    float gap;
    int numBlocks;
    int maxBlocks;
    float score_matrix[128][128];
};

#endif // ALIGN_H
