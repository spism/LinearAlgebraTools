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
     * vectors that map to the vector in the different coordinate system.  It calculates the reduced echelon form of the matrix.
     * @param augmentedMatrix is the input matrix with the solution vector included
     */
    public void findSolutionMatrix(double[][] augmentedMatrix)
    {
        int pivotCounter = 0;
        final double epsilon = 0.0000000001;
        for(int i = 0; i < augmentedMatrix[0].length; i++)
        {
            if(augmentedMatrix[0][i] == 0) // check for 0 in row 0, necessary because then the calculations are not possible and a solution can not be found
            {
                int nonZeroRow = -1;
                for(int j = 0; j < augmentedMatrix.length; j++)
                {
                    if(augmentedMatrix[j][i] != 0) nonZeroRow = j;
                }
                if(nonZeroRow > -1) interchangeRows(augmentedMatrix,pivotCounter,nonZeroRow);
                else continue; // no row that has nonzero value means this column does not need to be calculated, so increment column
            }
            for(int j = 0; j < augmentedMatrix.length; j++)
            {
                addRowMultiple(augmentedMatrix, -1 * (augmentedMatrix[j][i] / augmentedMatrix[0][i]),j,pivotCounter);
            }
            pivotCounter++;
        }
        for(int i = 0; i < augmentedMatrix.length; i++) // changes all pivots to 1
        {
            for(int j = 0; j < augmentedMatrix[0].length; j++)
            {
                if(augmentedMatrix[i][j] != 0)
                {
                    divideRow(augmentedMatrix,augmentedMatrix[i][j],i);
                    break;
                }
            }
        }
        pivotCounter = 0;
        for(int i = 0; i < augmentedMatrix[0].length; i++) // swaps rows so that first pivot is in first row, second in second, etc.
        {
            for(int j = 0; j < augmentedMatrix.length; j++)
            {
                if(augmentedMatrix[j][i] != 0)
                {
                    interchangeRows(augmentedMatrix,pivotCounter,j);
                    pivotCounter++;
                    break;
                }
            }
        }
        for(int i = 0; i < augmentedMatrix.length; i++) // calculates reduced echelon form of matrix
        {
            for(int j = 0; j < augmentedMatrix[0].length; j++)
            {
                if(augmentedMatrix[i][j] != 0)
                {
                    for(int k = i - 1; k >= 0; k--)
                    {
                        if(augmentedMatrix[k][j] != 0)
                        {
                            addRowMultiple(augmentedMatrix,-augmentedMatrix[k][j],k,i);
                        }
                    }
                    break;
                }
            }
        }
    }

    private void interchangeRows(double[][] matrix, int row1, int row2)
    {
        double[] tempRow = matrix[row2];
        matrix[row2] = matrix[row1];
        matrix[row1] = tempRow;
    }

    private void multiplyRow(double[][] matrix, double scalar, int row)
    {
        for(int i = 0; i < matrix[0].length; i++)
        {
            matrix[row][i] = matrix[row][i] * scalar;
        }
    }

    private void divideRow(double[][] matrix, double scalar, int row)
    {
        for(int i = 0; i < matrix[0].length; i++)
        {
            matrix[row][i] = matrix[row][i] / scalar;
        }
    }

    private void addRowMultiple(double[][] matrix, double scalar, int row, int rowToAdd)
    {
        for(int i = 0; i < matrix[0].length; i++)
        {
            matrix[row][i] += matrix[rowToAdd][i] * scalar;
        }
    }

    public double determinant2by2(double[][] matrix)
    {
        if(matrix.length > 2 || matrix[0].length > 2)
        {
            System.out.println("Matrix is not 2x2!");
            return 0;
        }
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    }

    public double determinant(double[][] matrix)
    {
        if(matrix.length == 2 && matrix[0].length == 2) return determinant2by2(matrix);

        double determinant = 0;
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

        int smallestZeroCount = matrix.length;
        int row = -1;
        int col = -1;
        for(int i = 0; i < matrix.length; i++)
        {
            if(rowZeroCounts[i] < smallestZeroCount)
            {
                smallestZeroCount = rowZeroCounts[i];
                row = i;
                col = -1;
            }
            if(colZeroCounts[i] < smallestZeroCount)
            {
                smallestZeroCount = colZeroCounts[i];
                row = -1;
                col = i;
            }
        }

        if(col > -1)
        {
            for(int i = 0; i < matrix.length; i++)
            {
                double[][] newMatrix = removeRow(matrix,i);
                newMatrix = removeCol(newMatrix,col);
                determinant += Math.pow(-1,i + 2 + col) * matrix[i][col] * determinant(newMatrix);
            }
        }
        else if(row > -1)
        {
            for(int i = 0; i < matrix[0].length; i++)
            {
                double[][] newMatrix = removeCol(matrix,i);
                newMatrix = removeRow(newMatrix,row);
                determinant += Math.pow(-1,row + i + 2) * matrix[row][i] * determinant(newMatrix);
            }
        }
        else
        {
            for(int i = 0; i < matrix.length; i++)
            {
                double[][] newMatrix = removeRow(matrix,i);
                newMatrix = removeCol(newMatrix,0);
                determinant += Math.pow(-1,i + 2) * matrix[i][0] * determinant(newMatrix);
            }
        }

        return determinant;
    }

    public double[][] removeRow(double[][] matrix, int row)
    {
        double[][] newMatrix = new double[matrix.length - 1][matrix[0].length];
        for(int i = 0; i < matrix.length; i++)
        {
            if(i < row)
            {
                for(int j = 0; j < matrix[0].length; j++)
                {
                    newMatrix[i][j] = matrix[i][j];
                }
            }
            if(i > row)
            {
                for(int j = 0; j < matrix[0].length; j++)
                {
                    newMatrix[i - 1][j] = matrix[i][j];
                }
            }
        }
        return newMatrix;
    }

    public double[][] removeCol(double[][] matrix, int col)
    {
        double[][] newMatrix = new double[matrix.length][matrix[0].length - 1];
        for(int i = 0; i < matrix[0].length; i++)
        {
            if(i < col)
            {
                for(int j = 0; j < matrix.length; j++)
                {
                    newMatrix[j][i] = matrix[i][j];
                }
            }
            if(i > col)
            {
                for(int j = 0; j < matrix.length; j++)
                {
                    newMatrix[j][i - 1] = matrix[i][j];
                }
            }
        }
        return newMatrix;
    }

    public double[] eigenvector(double[][] matrix, double eigenvalue)
    {

    }

    public double[] findEigenvalues(double[][] matrix) // matrix must be square as eigenvalues are only defined for square matrices
    {

    }
}
