public class functions {
    public double[] equzeros() {
        double[] zerosp1 = {0.0, 1.0};
        return zerosp1;
    }

    public double[] intnodes() {
        double[] intnodesp1 = {0.21132487, 0.78867513};
        return intnodesp1;
    }

    public double[] intweights() {
        double[] intweightsp1 = {0.5, 0.5};
        return intweightsp1;
    }

    public double philequdist(double xi, int node) {
        double ret;
        if(node==1) {
            ret = 1-xi; 
        }
        else if(node==2) {
            ret = xi; 
        }
        else {
            ret = 999999; 
            System.out.println("Unzulaessige Eingabe für Knoten");
        }
        return ret;
    }

    public double philequdistgrad(double xi, int node) {
        double ret;
        if(node==1) {
            ret = -1; 
        }
        else if(node==2) {
            ret = 1; 
        }
        else {
            ret = 999999; 
            System.out.println("Unzulaessige Eingabe für Knoten");
        }
        return ret;
    }

    // Funktion zur Multiplikation von Matrizen
    public double[][] matmult(double[][] mat_A, double[][] mat_B) {
        double[][] mat_AB = new double[mat_A.length][mat_B[0].length];
        for (int i=0; i<mat_A.length; i++) {
            for (int j=0; j<mat_B[0].length; j++) {
                for (int k=0; k<mat_A[0].length; k++ ) {
                    mat_AB[i][j] += mat_A[i][k]*mat_B[k][j];
                }
            }
        }
    return mat_AB;
    }

    // Funktion zur Multiplikation von Matrix und Skalar
    public double[][] matskal(double[][] mat, double skal) {
        for (int i=0; i<mat.length; i++) {
            for (int j=0; j<mat[0].length; j++) {
                mat[i][j] = mat[i][j] * skal;
            }
        }
    return mat;
    }

    // Funktion zur Addition von Matrizen, setzt gleiche Größe voraus
    public double[][] matplus(double[][] mat_A, double[][] mat_B) {
        double[][] mat_ApB = new double[mat_A.length][mat_A[0].length];
        for (int i=0; i<mat_A.length; i++) {
            for (int j=0; j<mat_A[0].length; j++) {
                mat_ApB[i][j] = mat_A[i][j] + mat_B[i][j];
            }
        }
    return mat_ApB;
    }

    public static void main(String[] args) {

    }
}
