/*************************************************************************
 *  Compilation:  javac LZWmod.java
 *  Execution:    java LZWmod - < input.txt > output.lzw  (compress input.txt 
 *                                                         into output.lzw)
 *  Execution:    java LZWmod + < output.lzw > input.rec  (expand output.lzw 
 *                                                         into input.rec)
 *  Dependencies: BinaryStdIn.java BinaryStdOut.java
 *
 *  Compress or expand binary input from standard input using LZW.
 *
 *
 *************************************************************************/

public class LZW {
    private static final int R = 256;        // alphabet size
    private static boolean flushIfFull = false;

    public static void compress() {
        CompressionCodeBookInterface codebook = 
            new DLBCodeBook(9, 16);

        // Writing the first bit of compressed file
        byte flushByte = 0;
        if(flushIfFull) flushByte = 1;
        BinaryStdOut.write(flushByte, 1);

        while (!BinaryStdIn.isEmpty()) {
            char c = BinaryStdIn.readChar();
            if(!codebook.advance(c)){ //found longest match
                int codeword = codebook.getCodeWord();
                BinaryStdOut.write(codeword, codebook.getCodewordWidth()); 
                codebook.add(flushIfFull);
                codebook.advance(c);
            }
        }
        int codeword = codebook.getCodeWord();
        BinaryStdOut.write(codeword, codebook.getCodewordWidth()); 

        BinaryStdOut.write(R, codebook.getCodewordWidth()); 
        BinaryStdOut.close();
    }


    public static void expand() {
        ExpansionCodeBookInterface codebook = new ArrayCodeBook(9, 16);

        // Read the first bit of compressed file
        flushIfFull = BinaryStdIn.readBoolean();

        int codeword = BinaryStdIn.readInt(codebook.getCodewordWidth(flushIfFull));
        String val = codebook.getString(codeword);

        while (true) {
            BinaryStdOut.write(val);
            codeword = BinaryStdIn.readInt(codebook.getCodewordWidth(flushIfFull));

            if (codeword == R) break;
            String s = codebook.getString(codeword);
            if (codebook.size() == codeword) s = val + val.charAt(0); // special case hack

            codebook.add(val + s.charAt(0), flushIfFull);
            val = s;

        }
        BinaryStdOut.close();
    }



    public static void main(String[] args) {
        if (args[0].equals("-")) {
          // Read the next arg to see if flushing the codebook will happen
          if (args[1].equalsIgnoreCase("n")) flushIfFull = false; // We won't flush the codebook
          else if (args[1].equalsIgnoreCase("r")) flushIfFull = true; // We will flush the codebook
          else throw new RuntimeException("Illegal command line argument");
          compress();
        }
        else if (args[0].equals("+")) expand();
        else throw new RuntimeException("Illegal command line argument");
    }

}
