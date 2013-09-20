using System;
using System.Drawing;
using System.Drawing.Imaging;
using System.Windows.Forms;
using System.Diagnostics;
using System.Collections.Generic;

class guiwin : Form
{
    Matrix originalImage;
    Matrix displayedImage;
    Matrix undoImage;

    OrientationField j;
    OrientationField jAve4;
    OrientationField jAve8;
    OrientationField jAve16;

    Matrix imageMask;


    Matrix Laplacian = new float [,] {
                              {  -1, -1, -1 },
                              {  -1,  8, -1 },
                              {  -1, -1, -1 }
                            };


    private static void clearBorder(Matrix image, int border, float pixel)
    {
        int height = image.Rows;
        int width = image.Columns;

        for (int r = 0; (r < height) && (r < border) ; r++)
            for (int c = 0; c < width; c++){
                image[r, c] = pixel;
                image[height - r - 1, c] = pixel;
            }

        for (int c = 0; (c < width) && (c < border) ; c++)
            for (int r = 0; (r < height) ; r++){
                image[r, c] = pixel;
                image[r, width - c - 1] = pixel;
            }
    }


    private Matrix filterGaussian(Matrix image, float variance, int filterRadius)
    {
        Matrix result;

        // naive full convolution
        // return ImageUtils.Convolve(image, ImageUtils.gaussianMask(4.0f, filterRadius) );

        // calculation of convolution with two 1-D filters
        Matrix kernelH = ImageUtils.gaussianRow(variance, filterRadius);

        result = ImageUtils.Convolve(image, kernelH);
        result = ImageUtils.Convolve(result, kernelH.Transpose);

        return result;
    }


        // toggle display between processed image and original
        // greyscale image.
    protected void toggle()
    {
        if (displayedImage != originalImage)
            displayImage(originalImage);
        else
            undo();
    }

    private void undo()
    {
        displayImage(undoImage);
    }

    private void displayImage(Matrix image)
    {
        undoImage = displayedImage;
        displayedImage = image;


        ClientSize = new Size(displayedImage.Columns, displayedImage.Rows);
        Invalidate();
    }

    private Matrix doImage(Matrix image)
    {
        originalImage = image;
        displayedImage = image;

        j = orientationField(image);
        jAve4 = filterField(j, 4);
        jAve8 = filterField(j, 8);
        jAve16 = filterField(j, 16);
    
        Matrix result = fingerprintEnhance(image);
        return result;
    }

    private void doFile(string filename)
    {
		Text = String.Format("Image:  {0}", filename);

        Matrix image = ImageUtils.GetGreyscale(filename);
        Matrix result = doImage(image);
        displayImage(result);
    }

	public guiwin(string filename)
	{
            // Set up drag and drop.
        AllowDrop = true;
        DragDrop += new
               System.Windows.Forms.DragEventHandler(dragDrop);
        DragEnter += new
               System.Windows.Forms.DragEventHandler(dragEnter);

        DoubleBuffered = true;

        doFile(filename);
	}

    


    double rotationDegrees = 15;
    double rotationIncrement = 15;

    protected override void OnKeyPress(KeyPressEventArgs e)
    {
        Matrix image = displayedImage;
        Matrix result = null;
        
        switch(e.KeyChar){
            case 'R':
                result = rotateImage(originalImage, rotationDegrees);
                rotationDegrees = (rotationDegrees + rotationIncrement) % 360;
                break;

            case 'a': result = autoScaleLevels(image); break;

            case 'f': result = fingerprintEnhance(image); break;
            case '1': result = experiment(image); break;                      
                      
            case 'm': imageMask = orientationMask(image); return;
            case 'p': result = orientationPrune(image); break;

            case 'r': lines.Clear(); result = restore(); break;

            case 'd': result = downscale(image, 2.5f); break;
            case 'D': result = downsample(image); break;
            case 'U': result = upscale(image); break;
            case 'i': result = invert(image); break;
            case 't': toggle(); return; // break;
            case 'T': result = threshold(image); break;
            case 's': save(image, "output.png"); return; // break;

            case 'g': result = filterGaussian(image, 1.0f, 8); break;

            case 'F': result = inRelief(image); break;


            case 'x': undo(); return; //  break;

            case 'q': quit(); return; // break;

            default:
                // no valid command typed.
                return;
        }
        displayImage(result);
    }



    protected void save(Matrix image, string filename)
    {
        ImageUtils.PutImage(filename, image);
        Console.WriteLine("image saved to file: {0}", filename);
    }


    protected Matrix autoScaleLevels(Matrix image)
    {
        Matrix r = image.Clone();
        ImageUtils.AutoScaleLevels(r);
        return r;
    }


    protected Matrix orientationMask(Matrix image)
    {
        Matrix result = image.Clone();

        image.Apply(delegate(int r, int c){
            double n = mag(jAve8, r, c);
            if ( n > 7.5){
                result[r, c] = 1;
            } else {
                result[r, c] = 0;
            }
           }
        );
        return result;
    }

    protected Matrix orientationPrune(Matrix image)
    {
        Matrix result = image.Clone();

        if (imageMask == null){
            Console.WriteLine("no image mask is defined");
            return result;
        }

        imageMask.Apply(delegate(int r, int c){
            double n = imageMask[r, c];
            if ( n > 0){
                result[r, c] = image[r, c];
            } else {
                result[r, c] = 255;
            }
           }
        );
        return result;
    }


    protected Matrix markupImage(Matrix image, OrientationField narrow, OrientationField wide)
    {
        Matrix result = image.Clone();

        image.Apply(delegate(int r, int c){
            double n = mag(narrow, r, c);
            double w = mag(wide, r, c);
            if ( (n > 50) &&  (n < 0.7 * w) ){
                result[r, c] = image[r, c];
            } else {
                result[r, c] = 0;
            }
           }
        );
        return result;
    }


    protected Matrix experiment(Matrix image)
    {
        int i;
        
//      doImage(displayedImage);

        for (i = 0; i < 10; i++){
            displayImage(fingerprintEnhance(displayedImage));
            Console.WriteLine("pass: {0}", i);
        }
        return displayedImage;
    }


    protected Matrix fingerprintEnhance(Matrix image)
    {
        Matrix result;

        const int filterRadius = 8;
        const float filterVariance = 4.0f;  // std deviation squared

//        const int filterRadius = 6;
//        const float filterVariance = 4.0f;  // std deviation squared


        result = filterGaussian(image, filterVariance, filterRadius);
        result = ImageUtils.Convolve(result, Laplacian);

        // We leave a few levels of grey.
        ImageUtils.BoundImage(result, -4.0f, -1.0f);

        ImageUtils.AutoScaleLevels(result);
        clearBorder(result, filterRadius + 1, 255); 

        // Uncomment the following three lines to get a pure b&w image.
        ImageUtils.InvertImage(result);
        ImageUtils.ThresholdImage(result, 1);
        ImageUtils.InvertImage(result);

        return result;
    }


    protected Matrix invert(Matrix image)
    {
        Matrix r = image.Clone();
        ImageUtils.InvertImage(r);
        return r;
    }

    protected Matrix threshold(Matrix image)
    {
        Console.WriteLine("threshold");
        Matrix r = image.Clone();
        ImageUtils.ThresholdImage(r, 128);
        return r;
    }




    // Interpolate a pixel value at floating coordinates r and c, using
    // the four nearest integral neighbors and the interpolating
    // polynomial in Lagrangian form.
    private float interpolate(Matrix m, double r, double c)
    {
        int rL = (int) r;
        int rU = rL + 1;
        int cL = (int) c;
        int cU = cL + 1;

        if ( (rL < 0) || ( cL < 0) )
            return 0;

        if ( (rL >= m.Rows) || (rU >= m.Rows) )
            return 0;
        
        if ( (cL >= m.Columns) || (cU >= m.Columns) )
            return 0;

        double v;

         v = (cU - c) * (rU - r) * m[rL, cL] +
             (r - rL) * (cU - c) * m[rU, cL] +
             (c - cL) * (rU - r) * m[rL, cU] +
             (c - cL) * (r - rL) * m[rU, cU];

        return (float) v;
    }


    /* From the specified center [cRow, cCol], go out in a direction
     * specified by degrees, as indicated by [row, col] along the
     * new axes.  Interpolate the result.
     */
    private float interpolate(Matrix m, int cRow, int cCol, double degrees, float row, float col)
    {
        double theta = degrees * Math.PI / 180;
        double sinTheta = Math.Sin(theta);
        double cosTheta = Math.Cos(theta);

        double dx = row * sinTheta + col * cosTheta;
        double dy = row * cosTheta - col * sinTheta;

        return interpolate(m, cRow + dy, cCol + dx);
    }



    private Matrix rotateAboutPoint(Matrix image, int x, int y, int radius, double degrees)
    {
        Matrix v = new Matrix(2 * radius + 1, 2 * radius + 1);
        v.Apply(delegate(int r, int c) {
                    v[r, c] = interpolate(image, y, x, degrees, r - radius, c - radius);
                }
        );
        return v;
    }



    private Matrix rotateImage(Matrix image, double degrees)
    {
        int cy = image.Rows / 2;
        int cx = image.Columns / 2;

        return rotateAboutPoint(image, cy, cx, Math.Max(cy, cx), degrees);
    }



    class guiLine {
        Point center;
        double ccwDegrees;

        public guiLine(Point c, double a, int len){
            center = c;
            ccwDegrees = a;
        }
                
        private void drawLine(Graphics g, double angle, Pen color)
        {
            int length = 60;
            g.SmoothingMode = System.Drawing.Drawing2D.SmoothingMode.AntiAlias;

            double radians = angle * Math.PI / 180;
            int dx = (int) (length * Math.Cos(radians)) / 2;
            int dy = (int) (length * -Math.Sin(radians)) / 2;

            Size delta = new Size(dx, dy);
            Point from = center - delta;
            Point to = center + delta;
            g.DrawLine(color, from, to);
        }

        public void render(Graphics g)
        {
            if (center.X < 0)
                return;
         
            drawLine(g, ccwDegrees, Pens.Red);
            drawLine(g, ccwDegrees + 90, Pens.Yellow);

        }
    }


// List of lines to render in orientation investigation.  (At the
// moment I only draw one, but such was not always he case!
    List<guiLine> lines = new List<guiLine>();


    private Matrix inRelief(Matrix image)
    {
        Matrix hDiff = new float[,] {
            {-1, 1, 0}
        };

        Matrix vDiff = hDiff.Transpose;


        Matrix gradH = new float[,] {
            {  0, -2, 0},
            {  0,  0, 0},
            {  0,  2, 0}
        };

        Matrix gradV = gradH.Transpose;

        Matrix gX = ImageUtils.Convolve(image, -hDiff);
        Matrix gY = ImageUtils.Convolve(image, -vDiff);

        Matrix r = (gX + gY) / 2;
//      ImageUtils.BoundImage(r, 0, 255);
        ImageUtils.AutoScaleLevels(r);

//      ImageUtils.PutImage("gx.png", r);
        return r;
    }


/* Experimental analysis of the neighborhood around
 * the point (x, y).
 */
    private void analyzePoint(int x, int y)
    {
        int guiLineLength = 60;
        Matrix image = displayedImage;

        lines.Clear();
        float orientation = localOrientationDegrees(image, x, y);
        lines.Add(new guiLine(new Point(x, y), orientation, guiLineLength));
        Invalidate();
    }


    // I do not overload the multiplication operator for Matrices
    // with this because that is best reserved for conventional
    // Matrix multiplication.  I could put this on the ^ operator,
    // perhaps.
    Matrix termwiseMultiply(Matrix lhs, Matrix rhs)
    {
        // assumes same size
        Matrix v = new Matrix(lhs.Rows, lhs.Columns);
        v.Apply(delegate(int r, int c) { v[r, c] = lhs[r, c] * rhs[r, c]; } );
        return v;
    }

    protected class OrientationField
    {
        public Matrix Re;
        public Matrix Im;
    }


    private OrientationField orientationField(Matrix image)
    {
        OrientationField j = new OrientationField();

        Matrix vDiff = new float[,] {
            { 0},
            {-1},
            { 1}
        };
        
        Matrix hDiff = new float[,] { {0, -1, 1} };
            
            // horizontal component of gradient
        Matrix gX = ImageUtils.Convolve(image, hDiff);

            // vertical component of gradient
        Matrix gY = ImageUtils.Convolve(image, vDiff);

        // The complex quantiy gX - i gY is the gradient matrix. (The
        // imaginary part is negative because the y axis increases
        // downward, yet we still measure angles positive ccw.) For
        // averaging directions, we create a complex number j which is
        // the square of g.  This effectively doubles the angles and
        // makes 180-degree differing angles the same.

        j.Im = -2 * termwiseMultiply(gX, gY);
        j.Re = termwiseMultiply(gX, gX) - termwiseMultiply(gY, gY);
        return j;
    }


    OrientationField filterField(OrientationField f, int radius)
    {
// why does this affect the speed?  very strange.

//      const float variance = 1000.0f;
      const float variance = 100.0f;
//        const float variance = 3.0f;
//        const float variance = 0.33f;
        OrientationField jAve = new OrientationField();
        jAve.Re = filterGaussian(f.Re, variance, radius);
        jAve.Im = filterGaussian(f.Im, variance, radius);
        return jAve;
    }

    double mag(OrientationField f, int r, int c)
    {
        double a = f.Re[r, c];
        double b = f.Im[r, c];
        return Math.Sqrt(a*a + b*b);
    }

    double radians(OrientationField f, int r, int c)
    {
        return Math.Atan2(f.Im[r, c], f.Re[r, c]) / 2;
    }

    double degrees(OrientationField f, int r, int c)
    {
        return radians(f, r, c) * 180.0 / Math.PI;
    }

    double report(string tag, OrientationField f, int x, int y)
    {
        // degrees is along the normal, so we add 90 degrees to get the
        // orientation.
        double deg = degrees(f, y, x) + 90;

        Console.WriteLine("{4}: ({0}, {1}): {2:0.00} degrees, abs={3:0.00}",
                          x, y, deg, mag(f, y, x), tag);
        return deg;
    }


    private float localOrientationDegrees(Matrix image, int x, int y)
    {
        double degrees = report(" jAve8", jAve8, x, y);
        return (float) degrees;
    }


    protected Matrix downsample(Matrix image)
    {
        Matrix smaller = new Matrix(image.Rows / 2, image.Columns / 2);
        smaller.Apply(delegate(int r, int c) { smaller[r, c] = image[2 * r, 2 * c]; } );
        return smaller;
    }

    protected Matrix downscale(Matrix img, float variance)
    {
        Console.WriteLine("downscale({0})", variance);

        // We apply a low-pass filter before downsampling.
        img = filterGaussian(img, variance, 8);
        Matrix smaller = downsample(img);
        return smaller;
    }

    protected Matrix upscale(Matrix img)
    {
        Matrix bigger = new Matrix(2 * img.Rows, 2 * img.Columns);

        bigger.Apply(delegate(int r, int c) { bigger[r, c] = img[r / 2, c / 2]; } );
        return bigger;
    }
    protected void quit()
    {
        Console.WriteLine("Bye!");
        Application.Exit();
    }


            // start over with the original image
    protected Matrix restore()
    {
        Console.WriteLine("start over with the original image");
        Matrix img = originalImage;
        return img;
    }

    protected override void OnMouseDown(MouseEventArgs e)
    {
        int x = e.X;
        int y = e.Y;

        if((ModifierKeys & Keys.Control) == Keys.Control){
            analyzePoint(x, y);
            return;
        } 

        if (e.Button == MouseButtons.Left){
            toggle();
        }
    }

	protected override void OnPaint( PaintEventArgs e )
	{
        Bitmap toDisplay;
            toDisplay = ImageUtils.PutImage(displayedImage);

        e.Graphics.DrawImage(
                             toDisplay,
                             new Rectangle(0, 0, toDisplay.Width, toDisplay.Height)
				     );

        foreach (guiLine line in lines){
            line.render(e.Graphics);
        }

	}

    private void dragEnter(object sender, System.Windows.Forms.DragEventArgs e)
    {
        if(e.Data.GetDataPresent(DataFormats.FileDrop))
            e.Effect = DragDropEffects.All;
        else
            e.Effect = DragDropEffects.None;
    }

    private void dragDrop(object sender, System.Windows.Forms.DragEventArgs e)
    {
        string[] s = (string[]) e.Data.GetData(DataFormats.FileDrop, false);
        int i;
        for(i = 0; i < s.Length; i++){
            doFile(s[i]);
        }
    }



}

class app
{
    [STAThread]
	public static void Main(string[] args)
	{
        string filename = "fingerx10.jpg";
        if (args.Length >= 1)
            filename = args[0];
		Application.Run(new guiwin(filename));
	}
}

