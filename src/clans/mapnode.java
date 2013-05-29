package clans;

public class mapnode {
    
    String level = "undefined";
    String info = "undefined";
    mapnode[] child = null;
    String leaf = "undefined";

    /**
     * see if this node contains a specific child
     * 
     * @param check
     * @return
     */
    public int haschild(mapnode check){
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
    public void add(mapnode addchild){
        if(child==null){
            child=new mapnode[1];
            child[0]=addchild;
        }else{
            int size = child.length;
            mapnode[] tmp=child;
            child=new mapnode[size+1];
            System.arraycopy(tmp,0,child,0,size);
            child[size]=addchild;
        }
    }   
}