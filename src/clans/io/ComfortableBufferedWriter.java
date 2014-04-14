package clans.io;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.Writer;

/**
 * A slightly more comfortable BufferedWriter to make saving runs with {@code ClusterData.saverun} easier to code.
 */
public class ComfortableBufferedWriter extends BufferedWriter {

	public ComfortableBufferedWriter(Writer out) {
		super(out);
		// TODO Auto-generated constructor stub
	}

	public ComfortableBufferedWriter(Writer out, int sz) {
		super(out, sz);
		// TODO Auto-generated constructor stub
	}
	
	public void print(Object in) throws IOException {
		write(in.toString());
	}
	
	public void println() throws IOException {
		write("\n");
	}
	
	public void println(Object in) throws IOException {
		write(in.toString() + "\n");
	}

}
