import java.util.Arrays;
import Jama.Matrix;

public class Main {

    // Elementsteifigkeitsmatrix
    static double[][] elematsteif(double length, double k) {
        functions getvalues = new functions();
        double[] intnodes = getvalues.intnodes();
        double[] intweights = getvalues.intweights();
        double philequdistgrad_i, philequdistgrad_ii;

        double[][] elesteifmat = {{0,0},{0,0}}; 
        double val_steif, val_mass;

        for(int j=0; j<=1; j++) {
            for(int i=0; i<=1; i++) {
                for(int ii=0; ii<=1; ii++) {
                    philequdistgrad_i = getvalues.philequdistgrad(intnodes[j], i+1);
                    philequdistgrad_ii = getvalues.philequdistgrad(intnodes[j], ii+1);

                    val_steif = k*(1/length)*philequdistgrad_i*philequdistgrad_ii*intweights[j];
                    elesteifmat[i][ii] = elesteifmat[i][ii]+val_steif;

                }
            }
        }
    return elesteifmat;
    }

    // Elementmassenmatrix
    static double[][] elematmass(double length, double rho, double cp) {
        functions getvalues = new functions();
        double[] intnodes = getvalues.intnodes();
        double[] intweights = getvalues.intweights();
        double philequdist_i, philequdist_ii;

        double[][] elemassmat = {{0,0},{0,0}};   
        double val_steif, val_mass;

        for(int j=0; j<=1; j++) {
            for(int i=0; i<=1; i++) {
                for(int ii=0; ii<=1; ii++) {
                    philequdist_i = getvalues.philequdist(intnodes[j], i+1);
                    philequdist_ii = getvalues.philequdist(intnodes[j], ii+1);

                    val_mass = rho*cp*length*philequdist_i*philequdist_ii*intweights[j];
                    elemassmat[i][ii] = elemassmat[i][ii]+val_mass;

                }
            }
        }
    return elemassmat;
    }

    // Assemblierung Systemmatrizen
    static double[][] assemble(double[][] elemat, int nrofelements) {
        int arraysize = nrofelements+1;
        double[][] sysmat = new double[arraysize][arraysize];
        int sysindex = 0;
        for(int i=0; i<=arraysize-1; i++) {
            for(int j=0; j<=arraysize-1; j++) {
                sysmat[i][j] = 0;
                        }
        }

        for(int ele=1; ele<=nrofelements; ele++) {
            sysmat[sysindex][sysindex] = sysmat[sysindex][sysindex]+elemat[0][0];
            sysmat[sysindex+1][sysindex] = sysmat[sysindex+1][sysindex]+elemat[1][0];
            sysmat[sysindex][sysindex+1] = sysmat[sysindex][sysindex+1]+elemat[0][1];
            sysmat[sysindex+1][sysindex+1] = sysmat[sysindex+1][sysindex+1]+elemat[1][1];
            sysindex++;
        }

    return sysmat;
    }

    // Zeitdiskretisierung
    static double[][] timeloop(double[][] mmat, double[][] kmat, double init, double dt, int ndt, double[] vector_rbs) {
        
        functions functions = new functions();        
        int len = mmat.length;
        //Copy array (of arrays...copyOf funktioniert daher nicht)
        double[][] mmatred = new double[len][];
            for(int i=0; i<len; i++) {
                mmatred[i] = Arrays.copyOf(mmat[i], len);
            }

        double gamma_l = vector_rbs[0];
        double gamma_r = vector_rbs[2];
        double g_l = vector_rbs[1];
        double g_r = vector_rbs[3];
        mmatred[0][0] -= gamma_l;
        mmatred[len-1][len-1] -= gamma_r;

        // Erstelle pmat und fülle mit ICs
        double[][] pmat = new double[len][1];
            for (int i=0; i<len; i++) {
                pmat[i][0] = init; 
            }

        double[][] sol_final = new double[len][1];
        //Timeloop
        for(int i=0; i<ndt; i++) {
            double[][] dtkmat = functions.matskal(kmat, -dt);
            double[][] mmat_dtkmat = functions.matplus(mmat, dtkmat);
            double[][] mmat_dtkmat_pmat = functions.matmult(mmat_dtkmat, pmat);
            double[][] rsvek = mmat_dtkmat_pmat;
            
            // Einbau RBs Lastvektor
            rsvek[0][0] -= g_l;
            rsvek[len-1][0] -= g_r;
            
            Matrix A = new Matrix(mmatred);
            Matrix B = new Matrix(rsvek);
            Matrix X = A.solve(B);

            double[][] solution = X.getArray();
                for (int ii = 0; ii < solution.length; ii++) {
                    pmat[ii][0] = solution[ii][0];
                    sol_final[ii][0] = solution[ii][0];
                 }         

        }

    return sol_final;
    }


    public static void main(String[] args) {
    
    // Materialdaten und Gesamtlänge
    double length = 0.1;
    double k = 2.5;
    double rho = 1000;
    double cp = 1000;
        
    // Randbedingungen
    double t_l = 25;    // Temperatur links 
    double h_l = 0.05;   // Widerstand links
    double t_r = -5;    // Temperatur rechts
    double h_r = 1;     // Widerstand rechts
    
    // Umrechnung Randbedingungen
    h_l = 1/h_l;
    h_r = 1/h_r;
    double gamma_l = -h_l;
    double gamma_r = -h_r;
    double g_l = -t_l*h_l;
    double g_r = -t_r*h_r;
    double[] vector_rbs = {gamma_l, g_l, gamma_r, g_r};

    // Anfangsbedingung
    double init = 0;
    double dt = 0.5;
    int ndt = 500;

    // Diskretisierung
    int nrofelements = 3;

    functions getvalues = new functions();
    double[] intnodes = getvalues.intnodes();
    double[] intweights = getvalues.intweights();
    double philequdist = getvalues.philequdist(0.35, 1);
    double philequdistgrad = getvalues.philequdistgrad(0.35, 1);
    
    double[][] elesteifmat = elematsteif(length/nrofelements, k);
    double[][] elemassmat = elematmass(length/nrofelements, rho, cp);


    double[][] sysmatsteif = assemble(elesteifmat, nrofelements);
    double[][] sysmatmass = assemble(elemassmat, nrofelements);

    System.out.println("Elementsteifigkeitsmatrix:");
    System.out.println(Arrays.deepToString(elesteifmat));
    System.out.println("Elementmassenmatrix:");
    System.out.println(Arrays.deepToString(elemassmat));

    if(nrofelements<=10) {
    System.out.println("Systemsteifigkeitsmatrix:");
    System.out.println(Arrays.deepToString(sysmatsteif));
    System.out.println("Systemmassenmatrix:");
    System.out.println(Arrays.deepToString(sysmatmass));
    }

    double[][] sol_final = timeloop(sysmatmass, sysmatsteif, init, dt, ndt, vector_rbs);
    System.out.println("Loesungswerte:");
    System.out.println(Arrays.deepToString(sol_final));    

    }

}
