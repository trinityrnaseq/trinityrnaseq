
public class BFLY_GLOBALS {

	public static int VERBOSE_LEVEL = 1;
	
	
	public static void debugMes(String mes, int verbosityLevel)
	{
		//TODO: use general logging that can be leveraged across all classes.

		if (verbosityLevel<=BFLY_GLOBALS.VERBOSE_LEVEL)
		{
			System.err.println(mes);
		}

	}

	
}
