#include <iostream>
#include <vector>
#include <string.h>

// upsampling test

const int IOBUFFSIZE = 1024;

struct sound_header {
    int ndim;		// number of dimensions
    int nchan;		// nummber of channels
    int len;    	// length of first dimension
    int fs;			// length of second dimension or sample rate
    int empty;		// length of third dimension
};

int main(int argc, char* argv[]){

    // determine impulse resopnse
    int L = 8;
    int D = 4;
    std::vector<float> h = {0.08, 0.25 , 0.64 , 0.95 , 0.95 , 0.64 , 0.25 , 0.08};
//    std::vector<float> x = {0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00 , 0.00};

    // read in file
    char* file_in = "../media/galway11_mono_downsampled2.bin";
    char* file_out = "../media/galway11_mono_downsampled.bin";

    struct sound_header header_x, header_y;

    FILE *fx, *fy;
    if (NULL == (fx = fopen(file_in, "rb"))) {
        printf("Error: Cannot find file for input.\n");
        return -1;
    }
    if (NULL == (fy = fopen(file_out, "wb"))) {
        printf("Error: Cannot open file for ouput.\n");
        return -1;
    }

    // read in header
    fread(&header_x, sizeof(header_x), 1, fx);

    if(header_x.nchan > 1) { printf("Error: function only takes 1 channel audio"); return -1; }
    std::cout << header_x.len << std::endl;
    std::cout << header_x.fs << std::endl;

    memcpy(&header_y, &header_x, sizeof(header_x));
    header_y.len = header_x.len/D+1;
    header_y.fs = header_x.fs/D;
    fwrite(&header_y, sizeof(header_x), 1, fy);

    // Processing
    float t; // Variable for accumulating convolution result
    int xlen, ylen = 0; // Indexes for input and output buffers
    int i;          // Index for input data buffer
    int l=1;        // Down sampling counter
    int k=0;        // Index for circular data buffer
    int n;          // Convolution loop index for filter coefficients and circular data buffer
    int test_sum = 0;
    float x[IOBUFFSIZE], y[IOBUFFSIZE];
    xlen = fread(x,sizeof(float),IOBUFFSIZE,fx); // Read in first chunk of input

    float *g = (float*)calloc(h.size(), sizeof(float));

    while(xlen>0) { // While there are samples to be processed, keep processing
        for(i=0; i<xlen; i++) { // Process each of the input samples
            k = (k+h.size()-1) % h.size();    // Update circular index of filter circular data buffer
            g[k] = x[i];          // Put each sample into the filter circular data buffer
            l = (l+D-1) % D;      // Update downsampling index
            if(l==0) {            // Down sampling condition: Compute convolution result when needed
                for(t=0.0, n=0; n<h.size(); n++) { // Convolution loop
                    t += h[n]*g[(n+k) % h.size()];   // Multiply and accumulate into local variable
                }                            // End convolution loop
                y[ylen] = t;           // Save convolution result into output buffer
                ylen++;                // Increment the index for the output buffer
                if(ylen==IOBUFFSIZE) { // If output buffer is full, then save it to output file
                    fwrite(y,sizeof(float),ylen,fy); // Write the output buffer
                    test_sum += ylen;
                    ylen = 0;               // Reset the index for the output buffer
                }
            } // End down sampling condition
        } // End loop over input samples
        xlen = fread(x,sizeof(float),IOBUFFSIZE,fx); // Read in next chunk of input samples
    }
    // Finish writing the last chunk of output samples
    if(ylen>0) {                       // If output buffer is full, then save it to output file
        fwrite(y,sizeof(float),ylen,fy); // Write the output buffer
        test_sum += ylen;
        ylen = 0;                        // Reset the index for the output buffer
    }
    std::cout << "\n" << test_sum << std::endl;
    fclose(fx);
    fclose(fy);
    free(g);
    std::cout << "success" <<std::endl;

    return 0;
}