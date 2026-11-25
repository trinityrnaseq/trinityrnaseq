// Source for ExitTester.jar. Trinity runs it to check that java starts (and can
// get the requested heap) before launching Butterfly: it exits with the given status.
//   usage: java -jar ExitTester.jar exitValue
// Rebuild with 'make jars' from the Trinity top-level directory.
public class ExitTester {

    public static void main(String[] args) {

        String usage = "\n\nusage: java -jar ExitTester.jar exitValue\n\n";

        if (args.length < 1) {
            System.err.print(usage);
            System.exit(1);
        }

        Integer exitValue = Integer.parseInt(args[0]);
        System.exit(exitValue.intValue());
    }
}
