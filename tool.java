public class tool
{
    public double[] scalarVector(double scalar, double[] vector) // assumes each vector is transposed
    {
        for(int i = 0; i < vector.length; i++) vector[i] = vector[i] * scalar;
        return vector;
    }

    public long[] scalarVector(long scalar, long[] vector)
    {
        for(int i = 0; i < vector.length; i++) vector[i] = vector[i] * scalar;
        return vector;
    }

    public double[] matrixVector(double[][] matrix, double[] vector)
    {
        double[] transformedVector = new double[matrix.length];
        for(int i = 0; i < matrix.length; i++)
        {
            for(int j = 0; j < matrix[i].length; j++)
            {
                transformedVector[i] += matrix[i][j] * vector[j];
            }
        }
        return transformedVector;
    }

    /*public double[] vectorVector(double[] vectorA, double[] vectorB)
    {
        double[] newVector = new double[vectorA.length];
        for(int i = 0; i < vectorA.length; i++)
        {
            newVector[i] = vectorA[i] * vectorB[i]; // no reason for this method to exist?
        }
        return newVector;
    }*/

    public double dotProduct(double[] vectorA, double[] vectorB)
    {
        double result = 0;
        for(int i = 0; i < vectorA.length; i++)
        {
            result += vectorA[i] * vectorB[i];
        }
        return result;
    }

    public double vectorNormal(double[] vector)
    {
        return Math.sqrt(dotProduct(vector,vector));
    }

    public double[][] matrixMatrix(double[][] matrixA, double[][] matrixB)
    {
        double[][] newMatrix = new double[matrixA.length][matrixB[0].length];
        for(int i = 0; i < matrixA.length; i++)
        {
            for(int j = 0; j < matrixB[0].length; j++)
            {
                double[] BVector = new double[matrixB.length];
                for(int k = 0; k < matrixB.length; k++)
                {
                    BVector[k] = matrixB[k][j];
                }
                newMatrix[i][j] = dotProduct(matrixA[i],BVector);
            }
        }
        return newMatrix;
    }

    public boolean isIdentity(double[][] matrix)
    {
        for(int i = 0; i < matrix.length; i++)
        {
            for(int j = 0; j < matrix[i].length; j++)
            {
                if((i != j && matrix[i][j] != 0) || matrix[j][j] != 1) return false;
            }
        }
        return true;
    }

    public double[][] inverse(double[][] matrix)
    {

    }

    /**
     * This method exists to find a solution matrix, which is calculated by creating an augmented matrix using the vector.  This is useful for finding possible
     * vectors that the vector will map to.
     * @param matrix is the input matrix
     * @param vector is the vector mapping in terms of the matrix
     * @return the matrix after row reduction to reduced echelon form
     */
    public double[][] findSolutionMatrix(double[][] matrix, double[] vector)
    {
        for(int i = 0; i < matrix[0].length; i++)
        {
            for(int j = 0; j < matrix.length; j++)
            {

            }
        }
    }

    private void interchange(double[][] matrix, int row1, int row2)
    {
        double[] tempRow = matrix[row2];
        matrix[row2] = matrix[row1];
        matrix[row1] = tempRow;
    }

    private void multiply(double[][] matrix, double scalar, int row)
    {
        for(int i = 0; i < matrix[0].length; i++)
        {
            matrix[row][i] = matrix[row][i] * scalar;
        }
    }

    private void addRowMultiple(double[][] matrix, double scalar, int row, int rowToAdd)
    {
        for(int i = 0; i < matrix[0].length; i++)
        {
            matrix[row][i] += matrix[rowToAdd][i] * scalar;
        }
    }

    public double determinant(double[][] matrix)
    {

    }

    public double[] eigenvector(double[][] matrix, double eigenvalue)
    {

    }

    public double[] findEigenvalues(double[][] matrix) // matrix must be square as eigenvalues are only defined for square matrices
    {
        int[] rowZeroCounts = new int[matrix.length];
        int[] colZeroCounts = new int[matrix.length];

        for(int i = 0; i < matrix.length; i++)
        {
            int colZeroCount = 0;
            int rowZeroCount = 0;
            for(int j = 0; j < matrix.length; j++)
            {
                if(matrix[i][j] == 0) rowZeroCount++;
                if(matrix[j][i] == 0) colZeroCount++;
            }
            rowZeroCounts[i] = rowZeroCount;
            colZeroCounts[i] = colZeroCount;
        }

        int smallestZeroCount = Integer.MAX_VALUE;
        int row;
        int col;
        for(int i = 0; i < matrix.length; i++)
        {
            if(rowZeroCounts[i] < colZeroCounts[i] && rowZeroCounts[i] < smallestZeroCount)
            {
                smallestZeroCount = rowZeroCounts[i];
            }
        }
    }
}
