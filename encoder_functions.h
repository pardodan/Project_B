#ifndef DSC_DNA_SEQ_ENCODER_FUNCTIONS_H
#define DSC_DNA_SEQ_ENCODER_FUNCTIONS_H
#include <vector>
#include <iostream>
using namespace std;

vector<vector<int>> split_read(vector<int>read , int k, int n);

unsigned long ** split_reed_solomon(vector<int>read ,int k , int n)



#endif //DSC_DNA_SEQ_ENCODER_FUNCTIONS_H
