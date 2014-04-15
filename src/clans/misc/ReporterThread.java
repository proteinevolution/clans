package clans.misc;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

/**
 * This class wraps a child thread and informs the parent about a completed child start and death.
 */
public class ReporterThread extends java.lang.Thread {

	private Object parent;
	private Thread child;

	private Object startLock;
	private Object runLock;
	
	private Method runOnCompletion;
	
	/**
	 * This class wraps a {@code child} thread and informs the {@code parent} about completed {@code child} start and
	 * death, i.e. the {@code child} is guaranteed to have isAlive() return true and false during the start completed
	 * and death completed notifications to the parent.
	 * <p>
	 * For informing about completed {@code child} thread start, the {@code start_lock} is notified. For informing about
	 * {@code child} thread death, the parameterless method @{code runOnCompletionMethodName} is invoked in the parent.
	 * <p>
	 * After starting the {@code child}, this thread waits for {@code run_lock.notify()} in the {@code child}, which
	 * should be called right before exiting.
	 * 
	 * @param parent
	 *            The parent object. This must contain a method with name {@code runOnCompletionMethodName}
	 * @param child
	 *            The child thread. This must call {@code run_lock.notify()} just before exiting.
	 * @param start_lock
	 *            The parent can use {@code start_lock.wait()} and will be signaled of a started child. Can be null if
	 *            the parent does not need this feature.
	 * @param run_lock
	 *            The child should be synchronized on this and then call {@code run_lock.notify()} just before returning
	 *            from its run() method.
	 * @param run_on_completion_method_name
	 *            Name of the method to be called in the parent when the {@code child} thread is not alive any more.
	 */
	public ReporterThread(Object parent, Thread child, Object start_lock, Object run_lock,
			String run_on_completion_method_name, String thread_name) {
		super();
		
		setName("ReporterThread");
		
		if (parent == null) {
			System.err.println("Error during construction of ReporterThread: "
					+ "received null as <parent> parameter");
			System.exit(101);
		}
		this.parent = parent;
		
		if (child == null) {
			System.err.println("Error during construction of ReporterThread: "
					+ "received null as <child> parameter");
			System.exit(101);
		}
		this.child = child;
		
		if (start_lock == null) {
			start_lock = new Object();
		}
		this.startLock = start_lock;
		this.runLock = run_lock;

		try {
			this.runOnCompletion = parent.getClass().getMethod(run_on_completion_method_name);
		
		} catch (NullPointerException e) {
			System.err.println("Error during construction of ReporterThread: "
					+ "received null as <run_on_completion_method_name> parameter");
			e.printStackTrace();
			System.exit(101);
			
		} catch (SecurityException e) {
			System.err.println("Error during construction of ReporterThread: "
					+ "Security manager does not allow getting method \"" + run_on_completion_method_name
					+ "\" from class \"" + parent.getClass() + "\n");
			System.err.println(e.getMessage());
			e.printStackTrace();
			System.exit(102);
			
		} catch (NoSuchMethodException e) {
			System.err.println("Error during construction of ReporterThread: "
					+ "Cannot instantiate ReporterThread as runOnCompletion method \"" + run_on_completion_method_name
					+ "\" was not found in class \"" + this.getClass() + "\"\n");
			System.err.println(e.getMessage());
			e.printStackTrace();
			System.exit(103);
		}
	}

	@Override
	public void run() {
		synchronized (runLock) {

			synchronized (startLock) {
				child.start();
				startLock.notify();
			}

			try {
				// this wait is canceled once the IterationsComputerThread.run method reaches its end
				runLock.wait();
			} catch (InterruptedException e) {
				System.err.println("fatal error: ReporterThread (child: " + child.getName()
						+ ") was interrupted during wait!\n");
				e.printStackTrace();
				System.exit(104);
			}
		}

		// make sure that the thread has finished all of its housekeeping an has actually died before continuing
		joinChild();
	}

	/**
	 * @return true if the {@code ReporterThread} is running and has an alive {@code child}.
	 */
	public boolean hasRunningIterationsComputerThread() {
		return isAlive() && child.isAlive();
	}

	/**
	 * Interrupts the child.
	 */
	public void stopChild() {
		child.interrupt();
	}

	/**
	 * Waits for the child thread to die and informs the parent about it.
	 */
	public void joinChild() {
		try {
			child.join(); // wait for the thread to finish

		} catch (InterruptedException e) {
			System.err.println("fatal error: ReporterThread (child: " + child.getName()
					+ ") was interrupted during joinChild!\n");
			e.printStackTrace();
			System.exit(105);
		}
		
		informParent();
	}
	
	/**
	 * Used to inform the parent once the child is not alive any more.
	 */
	private void informParent() {
		try {
			runOnCompletion.invoke(parent);

		} catch (IllegalArgumentException e) {
			System.err.println(e.getMessage());
			e.printStackTrace();
			System.exit(106);

		} catch (IllegalAccessException e) {
			System.err.println(e.getMessage());
			e.printStackTrace();
			System.exit(107);

		} catch (InvocationTargetException e) {
			System.err.println(e.getMessage());
			e.printStackTrace();
			System.exit(108);
		}
	}
}