using System;
using System.Drawing;
using System.Drawing.Imaging;
using System.Windows.Forms;
using System.Diagnostics;
using System.Runtime.InteropServices;

public class Matrix
{
    float[,] data;

    public Matrix(int rows, int cols)
    {
        data = new float[rows, cols];
    }

    public int Rows
    {
        get { return data.GetLength(0); }
    }
        
    public int Columns
    {
        get { return data.GetLength(1); }
    }


    // Construction of Matrix from array of floats.
    public Matrix(float[,] ar)
    {
        data = ar;
    }


    // automatic conversion from two-dimensional array of floats to
    // Matrix.
    public static implicit operator Matrix(float[,] d)
    {
        Matrix v = new Matrix(d);
        v.data = d;
        return v;
    }

    // automatic conversion from Matrix to two-dimensional array of
    // floats.
    public static implicit operator float[,](Matrix m)
    {
        return m.data;
    }


    public Matrix Transpose
    {
        get {
            Matrix v = new Matrix(Columns, Rows);
            v.Apply( delegate(int r, int c){ v.data[r, c] = this[c, r]; } );
            return v;
        }
    }


    // This lets us apply operations on a Matrix termwise, with a
    // minimum of tedious duplication.  It interoperates nicely with
    // anonymous methods and closures. See the implemention of
    // the arithmetic operators.
    
    public delegate void termfunc(int r, int c);
    public void Apply(termfunc func)
    {
        for (int r = 0; r < Rows; r++){
            for (int c = 0; c < Columns; c++){
                func(r, c);
            }
        }
    }


    // binary addition
    public static Matrix operator +(Matrix lhs, Matrix rhs) 
    {
        Matrix v = new Matrix(lhs.Rows, lhs.Columns);
        v.Apply( delegate(int r, int c) { v[r, c] = lhs[r, c] + rhs[r, c]; } );
        return v;
    }


    // binary subtraction
    public static Matrix operator -(Matrix lhs, Matrix rhs) 
    {
        Matrix v = new Matrix(lhs.Rows, lhs.Columns);
        v.Apply( delegate(int r, int c) { v[r, c] = lhs[r, c] - rhs[r, c]; } );
        return v;
    }

    // unary negation
    public static Matrix operator -(Matrix rhs) 
    {
        Matrix v = new Matrix(rhs.Rows, rhs.Columns);
        v.Apply(delegate(int r, int c) { v[r, c] = -rhs[r, c]; } );
        return v;
    }


    // scaling by a real number
    public static Matrix operator *(double lhs, Matrix rhs)
    {
        Matrix v = new Matrix(rhs.Rows, rhs.Columns);
        v.Apply(delegate(int r, int c) { v[r, c] = rhs[r, c] * (float) lhs; } );
        return v;
    }

    // multiplicative scaling is commutative
    public static Matrix operator *(Matrix lhs, double rhs)
    {
        return rhs * lhs;
    }

    // scaling by a real number via division
    public static Matrix operator /(Matrix lhs, double rhs)
    {
        Matrix v = new Matrix(lhs.Rows, lhs.Columns);
        v.Apply(delegate(int r, int c) { v[r, c] = lhs[r, c] / (float) rhs; } );
        return v;
    }

    public Matrix Clone()
    {
        Matrix copy = new Matrix(Rows, Columns);
        copy.Apply(delegate(int r, int c) { copy[r, c] = this[r, c]; } );
        return copy;
    }

    // array syntax for access
    public float this[int row, int col]
    {
        get {
            return data[row, col];
        }
        set {
            data[row, col] = value;
        }
    }

   public override string ToString()
    {
        string s = "";
        Apply(delegate(int r, int c){
            s += String.Format("{0} ", this[r, c]);
            if (c == Columns - 1)
                s += "\n";
         }
        );
  
      s += String.Format("\n");
        return s;
    }

}


static class ImageUtils
{

   public static string MatrixAsInts(Matrix m)
    {
        string s = "";
        m.Apply(delegate(int r, int c){
            s += String.Format("{0,5} ", (int) m[r, c]);
            if (c == m.Columns - 1)
                s += "\n";
         }
        );
  
      s += String.Format("\n");
        return s;
    }

   public static string MatrixAsBytes(Matrix m)
    {
        string s = "";
        m.Apply(delegate(int r, int c){
            s += String.Format("{0,3} ", (byte) m[r, c]);
            if (c == m.Columns - 1)
                s += "\n";
         }
        );
  
      s += String.Format("\n");
        return s;
    }

    public static void OffsetImage(Matrix image, float offset)
    {
        image.Apply(delegate(int r, int c){ image[r, c] += offset; } );
    }

    public static void InvertImage(Matrix image)
    {
        image.Apply(delegate(int r, int c){ image[r, c] = 255.0f - image[r, c]; } );
    }

    public static void BoundImage(Matrix image, float minLevel, float maxLevel)
    {
        image.Apply(delegate(int r, int c) {
            float pixel = image[r, c];
            if (pixel < minLevel)
                image[r, c] = minLevel;
            else if (pixel > maxLevel)
                image[r, c] = maxLevel;
        }
        );
    }


    public static void ThresholdImage(Matrix image, float whiteLevel)
    {
        image.Apply(delegate(int r, int c) {
            float pixel = image[r, c];
            if (pixel < whiteLevel)
                image[r, c] = 0;
            else
                image[r, c] = 255;
        }
        );
    }



    public static void AutoScaleLevels(Matrix image)
    {
        float max = float.MinValue, min = float.MaxValue;

        image.Apply(delegate(int r, int c) {
            float pixel = image[r, c];
            if (pixel > max)
                max = pixel;
            if (pixel < min)
                min = pixel;
        }
        );

        ScaleLevels(image, min, max);
    }

    public static void ScaleLevels(Matrix image, float blackLevel, float whiteLevel)
    {
        image.Apply(delegate(int r, int c) {
            float pixel = image[r, c];
            if (pixel >= whiteLevel)
                pixel = 255.0f;
            else if (pixel <= blackLevel)
                pixel = 0.0f;
            else
                pixel = (float) ( (pixel - blackLevel) * ( 256.0 / (whiteLevel - blackLevel)) )  ;
            image[r, c] = pixel;
        }
        );
    }


    public static Matrix gaussianRow(double variance, int radius)
    {
        int sideLength = 2 * radius + 1;
        Matrix m = new Matrix(1, sideLength);

        double sum = 0;
        for (int r = 0; r < sideLength; r++){
            for (int c = 0; c < sideLength; c++){
                int dx = (radius - c);
                int dy = (radius - r);
                sum += ( gaussian(dx, variance) * gaussian(dy, variance) );
            }
        }
        sum = Math.Sqrt(sum);

        for (int i = 0; i < sideLength; i++){
            int dx = (radius - i);
            m[0, i] = (float) (gaussian(dx, variance) / sum);
        }
        return m;
    }

    private static double gaussian(double x, double variance)
    {
        return Math.Exp( -(x * x) / (2 * variance) );
    }


    public static Matrix gaussianMask(double variance, int radius)
    {
        int sideLength = 2 * radius + 1;
        Matrix mask = new Matrix(sideLength, sideLength);

        double sum = 0.0f;
        for (int r = 0; r < sideLength; r++){
            for (int c = 0; c < sideLength; c++){
                int dx = (radius - c);
                int dy = (radius - r);
                double v = Math.Exp( -(dx*dx + dy*dy) / (2 * variance) );
                mask[r, c] = (float) v;
                sum += v;
            }
        }

        for (int r = 0; r < sideLength; r++){
            for (int c = 0; c < sideLength; c++){
                mask[r, c] /= (float) sum;
            }
        }
        return mask;
    }



    [DllImport("convolve.dll")]
            public static unsafe extern void convolve_gen(float *src, int height, int width,
                                                          float *coeff, int nRows, int nCols, float *dest);

    /* 
     * 2-D convolution.  This routine interfaces to a
     * C++ implementation in a DLL.
     */
    private static unsafe Matrix convolve_c(Matrix src, Matrix coeff)
    {
        int nRows = coeff.Rows;
        int nCols = coeff.Columns;

        int height = src.Rows;
        int width = src.Columns;

        Matrix dest = new Matrix(height, width);

        Stopwatch sw = new Stopwatch();
        sw.Start();

        fixed(float* dptr = (float[,]) dest){
            fixed(float* sptr = (float[,]) src){
                fixed (float *cptr = (float[,]) coeff){

                    // Call the C++ library version.
                    convolve_gen(sptr, height, width, cptr, nRows, nCols, dptr);

                }
            }
        }
        sw.Stop();
/*
        Console.WriteLine("convolve_c: image[{3}, {4}], mask[{1}, {2}], {0} milliseconds",
                          sw.Elapsed.TotalMilliseconds,
                          nRows, nCols,
                          src.Rows, src.Columns);
*/
        return dest;
    }

    // This is a version of convolve for masks that are column vectors.
    // The general convolve calls this for speed when appropriate.  No
    // similar version is needed for masks that are row vectors, because
    // the general convolve is already organized to exploit the data
    // horizontally.

    private static unsafe void convolve_column(float* src, int height, int width,
                                float* coeff, int nRows,
                                float* dest)
    {
        int ph = (nRows - 1) / 2;

        int h = height - nRows;
        int w = width - 1;

        float* dest_row = dest + (ph * width);

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


    /* 
     * 2-D convolution.  Self-contained version that does not call
     * anything outside this module.
     */
    private static unsafe float[,] convolve_net(Matrix src, Matrix coeff)
    {
        int height = src.Rows;
        int width = src.Columns;

        Matrix dest = new Matrix(height, width);

        Stopwatch sw = new Stopwatch();
        sw.Start();

        fixed(float* dptr = (float[,]) dest){
            fixed(float* sptr = (float[,]) src){
                fixed (float *cptr = (float[,]) coeff){

                    int nRows = coeff.Rows;
                    int nCols = coeff.Columns;

                    int ph = (nRows - 1) / 2;
                    int pc = (nCols - 1) / 2;

                    int stride = dest.Columns;

                    if (nCols == 1){
                        convolve_column(sptr, height, width, cptr, nRows, dptr);
                        return dest;  // verify: okay to return array that is pinned.
                    } 
                    float* row = sptr;
                    height -= nRows;
                    width -= nCols;
                    float* upper_left;

                    float* output_p = dptr + pc + (ph * stride);

                    for (int y = 0; y < height; y++, output_p += stride, row += stride){
                        float* outp = output_p;

                        upper_left = row;
                        for (int x = 0; x < width; x++, upper_left++){

                            float sum = 0.0f;  // faster than double here
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
                } // unfix coeff
            }  // unfix sptr
        } // unfix dptr

        sw.Stop();
        Console.WriteLine("convolve_net: image[{3}, {4}], mask[{1}, {2}], {0} milliseconds",
                          sw.Elapsed.TotalMilliseconds,
                          coeff.Rows, coeff.Columns,
                          src.Rows, src.Columns);

        return dest;
    }

    public static unsafe Matrix Convolve(Matrix src, Matrix coeff)
    {
        // We can call either the native C# version, or the C++ version
        // for a small speed improvement.
//      return convolve_net(src, coeff);
        return convolve_c(src, coeff);
    }



    public static void ExtractColorPlanes(Bitmap bitmap, out Matrix red, out Matrix green, out Matrix blue)
    {
        int rows = bitmap.Size.Height;
        int cols = bitmap.Size.Width;

        red =   new Matrix(rows, cols);
        green = new Matrix(rows, cols);
        blue =  new Matrix(rows, cols);

        // nyi: getting the bits

        int PixelSize = 0;

        switch (bitmap.PixelFormat)
        {
            case PixelFormat.Format24bppRgb:
                PixelSize = 3; break;
            case PixelFormat.Format32bppArgb:
                PixelSize = 4; break;
            case PixelFormat.Format32bppRgb:
                PixelSize = 4; break;
            case PixelFormat.Format8bppIndexed:
                PixelSize = 1; break;
            default:
                Console.WriteLine("Image format not supported");
                break;
        }

        unsafe {
            if (PixelSize == 1){

                // nyi: support for stupid indexed files simulating
                // greyscale
                
                throw new InvalidOperationException("Index files not supported.");

                
            } else if ((PixelSize == 3) || (PixelSize == 4) ){

                BitmapData bmd=bitmap.LockBits(
                                               new Rectangle(0, 0, bitmap.Width, bitmap.Height),
                                               ImageLockMode.ReadWrite,
                                               bitmap.PixelFormat
                                              );

                fixed(float* bluep = (float[,])blue){
                    fixed(float* greenp = (float[,])green){
                        fixed(float* redp = (float[,])red){

                            for(int y = 0; y < bmd.Height; y++) {
                                byte* row = (byte *) bmd.Scan0 + (y * bmd.Stride);

                                for(int x = 0; x < bmd.Width; x++){
                                    byte *pixel = &row[x * PixelSize];

                                    float alpha = 1.0f;
                                    if (PixelSize == 4)
                                        alpha = pixel[3] / 255.0f;

                                    bluep[y * cols + x] = alpha * pixel[0];
                                    greenp[y * cols + x] = alpha * pixel[1];
                                    redp[y * cols + x] = alpha * pixel[2];
                                }
                            }
                        }
                    }
                }

                bitmap.UnlockBits(bmd);

            } else if (PixelSize == 0){
                Console.WriteLine("PixelSize not determined");
            }
        }

    }

    public static Matrix GetGreyscale(Bitmap b)
    {
        Matrix red, green, blue;
        ExtractColorPlanes(b, out red, out green, out blue);
        return (red + green + blue) / 3;
    }

    public static Matrix GetGreyscale(string filename)
    {
        Bitmap b = new Bitmap(filename);
        return GetGreyscale(b);
    }




    public static unsafe Bitmap PutImage(Matrix red, Matrix green, Matrix blue)
    {
        int rows = red.Rows;
        int cols = red.Columns;

        if ( (blue.Rows != rows) || (green.Rows != rows) ||
             (blue.Columns != cols) || (green.Columns != cols) ){
                throw new InvalidOperationException("color plane sizes must all be the same");            
        }

        Bitmap bitmap = new Bitmap(cols, rows, PixelFormat.Format24bppRgb);

        BitmapData dest = bitmap.LockBits(
                                         new Rectangle(0, 0, bitmap.Width, bitmap.Height),
                                         ImageLockMode.WriteOnly,
                                         bitmap.PixelFormat
                                         );

        int stride = dest.Stride;

        fixed(float* bluep = (float[,]) blue){
            fixed(float* greenp = (float[,]) green){
                fixed(float* redp = (float[,]) red){
                    for (int r = 0; r < rows; r++)
                        for (int c = 0; c < cols; c++){
                            byte *p = ((byte *)dest.Scan0) + (r * stride) + 3*c;

                            p[0] = (byte) bluep[r * cols + c];
                            p[1] = (byte) greenp[r * cols + c];
                            p[2] = (byte) redp[r * cols + c];
                        }
                }
            }
        }

        bitmap.UnlockBits(dest);

        return bitmap;
    }


    public static Bitmap PutImage(Matrix grey)
    {
        return PutImage(grey, grey, grey);
    }

    public static void PutImage(string filename, Matrix grey)
    {
        Bitmap b = PutImage(grey);
        b.Save(filename);
    }



};  // class ImageUtils

