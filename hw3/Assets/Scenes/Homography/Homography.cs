using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using MathNet.Numerics.LinearAlgebra.Factorization;
using MathNet.Numerics.LinearAlgebra.Double;

public class Homography : MonoBehaviour {

    #region Public_Static_Methods

    public static double[,] CalcHomographyMatrix(double[,] s, double[,] d)
    {
        var x = CalcHomography(s, d);
        double[,] hGraphy = new double[3, 3];

        int row = 0;
        for (int i = 0; i < x.Length; i++)
        {
            if (i % 3 == 0 && i != 0)
                ++row;
            hGraphy[row, i % 3] = x[i];
        }
        return hGraphy;
    }

    public static double[,] CalcProjection(double[,] hGraphy, double[,] xy, bool log)
    {
        double[,] match = Homography.ApplyHomography(hGraphy, xy);
        if (log)
        {
            Debug.Log("...Homography Loading...");
            Debug.Log("(x,y) : " + xy[0, 0] + " , " + xy[1, 0] + " , " + xy[2, 0]);
            Debug.Log("(u,v) : " + match[0, 0] + " , " + match[1, 0] + " , " + match[2, 0]);
        }

        return match;
    }
    
    public static double[,] CalcInverseProjection(double[,] hGraphy, double[,] uv, bool log)
    {
        double[,] match = Homography.ApplyInverseHomography(hGraphy, uv);
        if (log)
        {
            Debug.Log("...Inverse Homography Loading...");
            Debug.Log("(u,v) : " + uv[0, 0] + " , " + uv[1, 0] + " , " + uv[2, 0]);
            Debug.Log("(x,y) : " + match[0, 0] + " , " + match[1, 0] + " , " + match[2, 0]);
        }

        return match;
    }

    #endregion

    #region Private_Static_Methods

    private static double[] CalcHomography(double[,] s, double[,] d)
    {
        int N = s.GetLength(0);
        int M = 9;

        double[,] Arr = new double[2 * N, M];

        
        int point = 0;
        for (int i = 0; i < 2 * N; i++)
        {
            if (i % 2 == 0)
            {
                Arr[i, 0] = -1 * s[point, 0];
                Arr[i, 1] = -1 * s[point, 1];
                Arr[i, 2] = -1;
                Arr[i, 3] = 0;
                Arr[i, 4] = 0;
                Arr[i, 5] = 0;
                Arr[i, 6] = s[point, 0] * d[point, 0];
                Arr[i, 7] = d[point, 0] * s[point, 1];
                Arr[i, 8] = d[point, 0];
            }
            else
            {
                Arr[i, 0] = 0;
                Arr[i, 1] = 0;
                Arr[i, 2] = 0;
                Arr[i, 3] = -1 * s[point, 0];
                Arr[i, 4] = -1 * s[point, 1];
                Arr[i, 5] = -1;
                Arr[i, 6] = s[point, 0] * d[point, 1];
                Arr[i, 7] = d[point, 1] * s[point, 1];
                Arr[i, 8] = d[point, 1];
                ++point;
            }
        }

        var mat = DenseMatrix.OfArray(Arr);
        var svd = mat.Svd(true);
        
        return svd.VT.Row(svd.VT.RowCount - 1).ToArray();

    }

    private static double[,] ApplyInverseHomography(double[,] hGraphy, double[,] a)
    {
        var m = DenseMatrix.OfArray(hGraphy).Inverse();
        return ApplyHomography(m.ToArray(), a);
    }

    private static double[,] ApplyHomography(double[,] hGraphy, double[,] a)
    {
        double[,] temp = MultiplyMatrices(hGraphy, a);
        double[,] res = new double[,] { { temp[0, 0] / temp[2, 0] }, { temp[1, 0] / temp[2, 0] }, { 1 } };

        return res;
    }

    public static double[,] MultiplyMatrices(double[,] a, double[,] b)
    {
        int m = a.GetLength(0);
        int n = b.GetLength(1);

        double[,] res = new double[m, n];

        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                res[i, j] = 0;
                for (int k = 0; k < a.GetLength(1); k++)
                {
                    res[i, j] += a[i, k] * b[k, j];
                }
            }
        }
        return res;
    }

    #endregion
}
