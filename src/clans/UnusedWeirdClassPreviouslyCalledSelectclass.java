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
public class UnusedWeirdClassPreviouslyCalledSelectclass {
    
    /** Creates a new instance of selectclass */
    public UnusedWeirdClassPreviouslyCalledSelectclass() {
    }
    
    public UnusedWeirdClassPreviouslyCalledSelectclass(int[] selectednames){
        this.selectednames=selectednames;
    }

    public UnusedWeirdClassPreviouslyCalledSelectclass(int[] selectednames, java.awt.Color color){
        this.selectednames=selectednames;
        this.color=color;
    }
    
    public UnusedWeirdClassPreviouslyCalledSelectclass(int[] selectednames, java.awt.Color color, String name){
        this.selectednames=selectednames;
        this.color=color;
        this.name=name;
    }
    
    String name="no name set";
    int[] selectednames=new int[0];
    java.awt.Color color=new java.awt.Color(0.5f,0.5f,0.5f);
}
