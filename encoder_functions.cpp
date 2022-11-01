#include "encoder_functions.h"

vector<vector<int>> split_read(vector<int>read ,int k , int n)
{
    vector<vector<int>> vec;


    int num_of_words = k/n + 1;

    // Inserting elements into vector
    for (int i = 0; i < num_of_words ; i++) {

        // Vector to store "column" elements
        vector<int> v1;

        for (int j = 0; j < n; j++) {
            if((n*i + j)<k) {
                v1.push_back(read[(n * i) + j]);
            }
            //padding with zeroes the last vector
            else{
                v1.push_back(0);
            }

        }

        // Pushing back above 1D vector
        // to create the 2D vector
        vec.push_back(v1);
    }

    cout << "The 2-D vector is:" <<endl;

    // Displaying the 2D vector
    for (int i = 0; i < vec.size(); i++) {
        for (int j = 0; j < vec[i].size(); j++)
            cout << vec[i][j] << " ";
        cout << endl;
    }
    return vec;
}

unsigned long ** split_reed_solomon(vector<int>read ,int k , int n)
{
//    RS_WORD *originalData = (RS_WORD*)malloc(sizeof(RS_WORD) * (rs_k + rs_nsym));
//    vector<vector<int>> vec;


    int num_of_words = k/n + 1;

//    // Inserting elements into vector
//    for (int i = 0; i < num_of_words ; i++) {
//
//        // Vector to store "column" elements
//        vector<int> v1;
//
//        for (int j = 0; j < n; j++) {
//            if((n*i + j)<k) {
//                v1.push_back(read[(n * i) + j]);
//            }
//                //padding with zeroes the last vector
//            else{
//                v1.push_back(0);
//            }
//
//        }
//
//        // Pushing back above 1D vector
//        // to create the 2D vector
//        vec.push_back(v1);
//    }
//
//    cout << "The 2-D vector is:" <<endl;
//
//    // Displaying the 2D vector
//    for (int i = 0; i < vec.size(); i++) {
//        for (int j = 0; j < vec[i].size(); j++)
//            cout << vec[i][j] << " ";
//        cout << endl;
//    }

    unsigned long** a = new unsigned long*[num_of_words];
    for(int i = 0; i < num_of_words; ++i) {
        a[i] = new unsigned long[n];
        for (int j = 0; j < n; j++) {
            if((n*i + j)<k) {

                a[i][j] = read[(n * i) + j];
                //v1.push_back(read[(n * i) + j]);
            }
                //padding with zeroes the last vector
            else{
                a[i][j] = 0;
                //v1.push_back(0);
            }

        }

    }
    return a;
}
