package clans;

public class MicroarrayMapnode {
    
    String level = "undefined";
    String info = "undefined";
    MicroarrayMapnode[] child = null;
    String leaf = "undefined";

    /**
     * see if this node contains a specific child
     * 
     * @param check
     * @return
     */
    public int haschild(MicroarrayMapnode check){
        if(child==null){
            return -1;
        }

        for (int i = child.length; --i >= 0;) {
            if(check==child[i]){
                return i;
            }
        }
        return -1;
    }
    
    /**
     * 
     * @param addchild
     */
    public void add(MicroarrayMapnode addchild){
        if(child==null){
            child=new MicroarrayMapnode[1];
            child[0]=addchild;
        }else{
            int size = child.length;
            MicroarrayMapnode[] tmp=child;
            child=new MicroarrayMapnode[size+1];
            System.arraycopy(tmp,0,child,0,size);
            child[size]=addchild;
        }
    }   
}