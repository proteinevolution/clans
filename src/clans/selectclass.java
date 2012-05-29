/*
 * selectclass.java
 *
 * Created on February 5, 2004, 10:18 AM
 */
package clans;
/**
 *
 * @author  tancred
 */
public class selectclass {
    
    /** Creates a new instance of selectclass */
    public selectclass() {
    }
    
    public selectclass(int[] selectednames){
        this.selectednames=selectednames;
    }

    public selectclass(int[] selectednames, java.awt.Color color){
        this.selectednames=selectednames;
        this.color=color;
    }
    
    public selectclass(int[] selectednames, java.awt.Color color, String name){
        this.selectednames=selectednames;
        this.color=color;
        this.name=name;
    }
    
    String name="no name set";
    int[] selectednames=new int[0];
    java.awt.Color color=new java.awt.Color(0.5f,0.5f,0.5f);
}
