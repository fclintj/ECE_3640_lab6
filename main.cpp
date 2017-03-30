#include <iostream>
#include <vector>
#include <string.h>
#include <cmath>

// upsampling test

const int IOBUFFSIZE = 1024;

typedef struct {
    int ndim;		// number of dimensions
    int nchan;		// nummber of channels
    int len;    	// length of first dimension
    int fs;			// length of second dimension or sample rate
    int empty;		// length of third dimension
} sound_header;

typedef struct {
    int ndim;		// number of dimensions
    int nchan;		// nummber of channels
    int len;    	// length of first dimension
    int empty1;			// length of second dimension or sample rate
    int empty2;		// length of third dimension
} transfer_function;

int main(int argc, char* argv[]){

    // determine impulse resopnse






    // read in file
    char* file_in = argv[1];
    char* file_out = argv[2];
    char* file_h = argv[3];
    int U = atoi(argv[4]);
    int D = atoi(argv[5]);

    int M;  // size of buffer


    sound_header header_x, header_y;
    transfer_function header_h;

    FILE *fx, *fy, *fh;
    if (NULL == (fx = fopen(file_in, "rb"))) {
        printf("Error: Cannot find file for input.\n");
        return -1;
    }
    if (NULL == (fh = fopen(file_h, "rb"))) {
        printf("Error: Cannot open transfer function file for processing.\n");
        return -1;
    }
    if (NULL == (fy = fopen(file_out, "wb"))) {
        printf("Error: Cannot open file for ouput.\n");
        return -1;
    }


    // read in header
    fread(&header_x, sizeof(header_x), 1, fx);
    fread(&header_h, sizeof(header_h), 1, fh);

    M = header_h.len;
    float *h = new float[M];
    fread(h,sizeof(float),M,fh);
    M = M-1;
    if(header_x.nchan > 1) { printf("Error: function only takes 1 channel audio"); return -1; }

    memcpy(&header_y, &header_x, sizeof(header_x));

    header_y.len = round(header_x.len*double(U)/D)+1;
    header_y.fs = header_x.fs*U/D;

    fwrite(&header_y, sizeof(header_x), 1, fy);

    float t;          // Variable for accumulating convolution result
    int xlen, ylen=0; // Indexes for input and output buffers
    int i;            // Index for input data buffer
    int j;            // Index for up sampling loop
    int k=0;          // Index for circular data buffer
    int m;            // Convolution loop indx for filter coefficients
    int n;            // Convolution loop index for circular data buffer
    int l=1;          // Down sampling counter

    int test_sum = 0;
    float x[IOBUFFSIZE], y[IOBUFFSIZE];
    float *g = (float*)calloc(M, sizeof(float));

    xlen = fread(x,sizeof(float),IOBUFFSIZE,fx);// Read in first chunk of input samples
    while(xlen>0) {                             // While there are samples to be processed, keep processing
        for(i=0; i<xlen; i++) {                 // Process each of the input samples
            k = (k+M-1) % M;      // Update circular index of filter circular data buffer
            g[k] = x[i];                        // Put each sample into the filter circular data buffer
            l = (l+D-1) % D;
            if(l==0) {                          // Down sampling condition: Compute convolution result when needed
                for (j = 0; j < U; j++) {       // Loop over the up sampled outputs
                    for (t = 0.0, m = 0, n = 0; n < M; n++, m += U) { // Convolution loop
                        t += h[m + j] * g[(n + k) % M];            // Multiply and accumulate into local variable
                    }                           // End convolution loop
                    y[ylen] = t;                // Save result into output buffer
                    ylen++;                     // Increment the index for the output buffer
                    if (ylen == IOBUFFSIZE) {   // If output buffer is full, then save it to output file
                        fwrite(y, sizeof(float), ylen, fy); // Write the output buffer
                        test_sum += ylen;
                        ylen = 0;               // Reset the index for the output buffer
                    }
                }                               // End loop over up sampled outputs
            }
        }                                       // End loop over input samples
        xlen = fread(x,sizeof(float),IOBUFFSIZE,fx); // Read in next chunk of input samples
    }
    // Finish writing the last chunk of output samples
    if(ylen>0) {                                // If output buffer is full, then save it to output file
        fwrite(y,sizeof(float),ylen,fy);        // Write the output buffer
        test_sum += ylen;
        ylen = 0;                               // Reset the index for the output buffer
    }
    fclose(fx);
    fclose(fy);
    free(g);
    std::cout << "header: " << header_y.len << std::endl;
    std::cout << "actual written: " << test_sum << std::endl;
    return 0;
}