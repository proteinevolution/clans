package clans.misc;

/**
 * Operating system (OS)-specific code.
 * <p>
 * This can be used to adjust for Java behaving differently on different OS, e.g. in certain File operations.
 * <p>
 * This code is from: http://stackoverflow.com/a/228499/454402
 */
public final class OsUtils {
	private static String OS = null;

	public static String getOsName() {
		if (OS == null) {
			OS = System.getProperty("os.name");
		}
		return OS;
	}

	public static boolean isWindows() {
		return getOsName().startsWith("Windows");
	}
	
	public static boolean isLinux() {
		return getOsName().startsWith("Linux");
	}
}