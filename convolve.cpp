extern "C"
{


// This is a version of convolve for masks that are column vectors.
// The general convolve calls this for speed when appropriate.  No
// similar version is needed for masks that are row vectors, because
// the general convolve is already organized to exploit the data
// horizontally.

void convolve_column(float *src, const int height, const int width,
                                        float *const coeff, const int nRows,
                                        float *const dest)
{
    int ph = (nRows - 1) / 2;

    const int h = height - nRows;
    const int w = width - 1;
    
    float *dest_row = dest + (ph * width);

    for (int y = 0; y < h; y++, dest_row += width, src += width){
        float* outp = dest_row;
        float* src_top = src;
        for (int x = 0; x < w; x++, src_top++){

            float sum = 0.0f;
            float *inp;
            int v;

            for (v = 0, inp = src_top; v < nRows; v++, inp += width){
                sum += ( *inp * coeff[v] );
            }
            // record the result
            *outp++ = sum;
        }
    }

}


__declspec(dllexport) void convolve_gen(float *src, const int height, const int width,
                                        float *coeff, const int nRows, const int nCols,
                                        float *dest)
{
    if (nCols == 1){
        return convolve_column(src, height, width, coeff, nRows, dest);
    }

    int ph = (nRows - 1) / 2;
    int pc = (nCols - 1) / 2;

    const int stride = width;

    const int h = height - nRows;
    const int w = width - nCols;
    
    float *cptr = coeff;
    float* row = src;
    float* output_p = dest + pc + (ph * width);

    for (int y = 0; y < h; y++, output_p += stride, row += stride){
        float* outp = output_p;
        float* upper_left = row;
        for (int x = 0; x < w; x++, upper_left++){

            float sum = 0.0f;  // faster than double here, at least in C#

            float* coeff_data = cptr;
            float* inp;
            int v;

            for (v = 0, inp = upper_left; v < nRows; v++, inp += stride){
                float* row_data = inp;
                for (int h = 0; h < nCols; h++){
                    sum += ( *row_data++ * *coeff_data++ );
                }
            }
            // record the result
            *outp++ = sum;
        }
    }
}

} // extern "C"
