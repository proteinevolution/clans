/*
 * polygons.java
 *
 * Created on June 8, 2005, 5:02 PM
 */
package clans;
/**
 *
 * @author  tancred
 */
public class makepolygons {
    
    /** Creates a new instance of polygons */
    public makepolygons() {
    }
    
    static int polygonnum=5;
    static float[][] square={{-1,1},{1,1},{1,-1},{-1,-1}};
    static float[][] diamond={{0,1.4f},{1.4f,0},{0,-1.4f},{-1.4f,0}};
    static float[][] downtriangle={{-1.5f,0.87f},{1.5f,0.87f},{0,-1.72f}};
    static float[][] uptriangle={{-1.5f,-0.87f},{1.5f,-0.87f},{0,1.72f}};
    static float[][] fourstar={{0,1.5f},{0.4f,0.4f},{1.5f,0},{0.4f,-0.4f},{0,-1.5f},{-0.4f,-0.4f},{-1.5f,0},{-0.4f,0.4f}};
    static float[][] fourstar90={{0,0.55f},{1,1},{0.55f,0},{1,-1},{0,-0.55f},{-1,-1},{-0.55f,0},{-1,1}};
    static float[][] sixstar={{0,1.15f},{0.33f,0.59f},{1,0.59f},{0.66f,0},{1,-0.59f},{0.33f,-0.59f},{0,-1},{-0.33f,-0.59f},{-1,-0.59f},{-0.66f,0},{-1,0.59f},{-0.33f,0.59f}};
    static float[][] upfivestar={{0,1.5f},{0.3f,0.4f},{1.4f,0.46f},{0.48f,-0.15f},{0.88f,-1.2f},{0,-0.5f},{-0.88f,-1.2f},{-0.48f,-0.15f},{-1.4f,0.46f},{-0.3f,0.46f}};
    static float[][] downfivestar={{0,-1.5f},{0.3f,-0.4f},{1.4f,-0.46f},{0.48f,0.15f},{0.88f,1.2f},{0,0.5f},{-0.88f,1.2f},{-0.48f,0.15f},{-1.4f,-0.46f},{-0.3f,-0.46f}};
    static float[][] plus={{0.4f,1.2f},{0.4f,0.4f},{1.2f,0.4f},{1.2f,-0.4f},{0.4f,-0.4f},{0.4f,-1.2f},{-0.4f,-1.2f},{-0.4f,-0.4f},{-1.2f,-0.4f},{-1.2f,0.4f},{-0.4f,0.4f},{-0.4f,1.2f}};
    static float[][] plus45={{0,0.57f},{0.57f,1.135f},{1.13f,0.57f},{0.57f,0},{1.13f,-0.57f},{0.57f,-1.13f},{0,-0.57f},{-0.57f,-1.13f},{-1.13f,-0.57f},{-0.57f,0},{-1.13f,0.57f},{-0.57f,1.13f}};
    static float[][] upY={{0,0.58f},{0.87f,1.08f},{1.37f,0.21f},{0.5f,-0.29f},{0.5f,-1.29f},{-0.5f,-1.29f},{-0.5f,-0.29f},{-1.37f,0.21f},{-0.87f,1.08f}};
    static float[][] downY={{0,-0.58f},{0.87f,-1.08f},{1.37f,-0.21f},{0.5f,0.29f},{0.5f,1.29f},{-0.5f,1.29f},{-0.5f,0.29f},{-1.37f,-0.21f},{-0.87f,-1.08f}};
    
    
    static java.util.ArrayList <int[][]> get(int size){
        //get all polygons defined in this method scaled to size
        java.util.ArrayList <int[][]>retvec=new java.util.ArrayList<int[][]>();
        //float mysize=((float)size)/2;
        retvec.add(new int[0][0]);//a polygon of type 0 should be an oval (default)
        retvec.add(getpol(square,size));
        retvec.add(getpol(diamond,size));
        retvec.add(getpol(uptriangle,size));
        retvec.add(getpol(downtriangle,size));
        retvec.add(getpol(fourstar,size));
        retvec.add(getpol(fourstar90,size));
        retvec.add(getpol(sixstar,size));
        retvec.add(getpol(upfivestar,size));
        retvec.add(getpol(downfivestar,size));
        retvec.add(getpol(plus,size));
        retvec.add(getpol(plus45,size));
        retvec.add(getpol(upY,size));
        retvec.add(getpol(downY,size));
        return retvec;
    }//end get
    
    static int[][] get(int number, int size){
        //get the array of int describing a specific polygon of a given size
        //System.out.println("in makepolygons for "+number+":"+size);
        if(number==0){
            return new int[0][0];//a polygon of type 0 should be an oval (default)
        }else if(number==1){
            return getpol(square,size);
        }else if(number==2){
            return getpol(diamond,size);
        }else if(number==3){
            return getpol(uptriangle,size);
        }else if(number==4){
            return getpol(downtriangle,size);
        }else if(number==5){
            return getpol(fourstar,size);
        }else if(number==6){
            return getpol(fourstar90,size);
        }else if(number==7){
            return getpol(sixstar,size);
        }else if(number==8){
            return getpol(upfivestar,size);
        }else if(number==9){
            return getpol(downfivestar,size);
        }else if(number==10){
            return getpol(plus,size);
        }else if(number==11){
            return getpol(plus45,size);
        }else if(number==12){
            return getpol(upY,size);
        }else if(number==13){
            return getpol(downY,size);
        }else{
            System.err.println("unknown polygon number "+number);
            return new int[0][0];
        }
    }//end get
    
    static int[][] getpol(float[][] inarr, float size){
        int elements=java.lang.reflect.Array.getLength(inarr);
        int[][] retarr=new int[3][elements];
        for(int i=0;i<elements;i++){
            retarr[0][i]=(int)(inarr[i][0]*size/2);
            retarr[1][i]=(int)(inarr[i][1]*size/2);
            retarr[2][i]=elements;
        }//end for i
        return retarr;
    }//end getpol
    
}
