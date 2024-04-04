/**
 * An implementation of ExpansionCodeBookInterface using an array.
 */

public class ArrayCodeBook implements ExpansionCodeBookInterface {
    private static final int R = 256;        // alphabet size
    private String[] codebook;
    private int W;       // current codeword width
    private int minW;    // minimum codeword width
    private int maxW;    // maximum codeword width
    private int L;       // maximum number of codewords with 
                         // current codeword width (L = 2^W)
    private int code;    // next available codeword value
  
    public ArrayCodeBook(int minW, int maxW){
        this.maxW = maxW;
        this.minW = minW;
        initialize();    
    }
    public int size(){
        return code;
    }


    public int getCodewordWidth(boolean flushIfFull){ 
      if(code == L) {
        // Maximum number of codewords has been reached
        if(W < maxW) {
          // We can add more codewords to the codebook; increase the size
          return W+1;
        }else if(W == maxW) {
          // We've hit the maximum codebook size; should we flush?
          if(flushIfFull) {
            // Yes
            return minW;
          } else {
            // No
            return maxW;
          }
        }
      }
      return W;
    }
    
    private void initialize(){
        codebook = new String[1 << maxW];
        W = minW;
        L = 1<<W;
        code = 0;
        // initialize symbol table with all 1-character strings
        for (int i = 0; i < R; i++)
            add("" + (char) i, false);
        add("", false); //R is codeword for EOF
    }

    public void add(String str, boolean flushIfFull) {
        boolean haveRoom = false;
        if(code < L){
          // Codebook still not full
            haveRoom = true;
        } else {
          // Codebook is full; what do we do?
          if (W < maxW) {
            // Haven't hit the max codebook expansion size
            W++;
            L = L*2;
            haveRoom = true;
          } else if (W == maxW) {
            // We have hit the max; do we flush?
            if (flushIfFull) {
              // Yes
              initialize();
              haveRoom = true;
            }
          }
        }

        if(haveRoom){
            codebook[code] = str;
            code++;
        }
    }

    public String getString(int codeword) {
        return codebook[codeword];
    }
    
}
