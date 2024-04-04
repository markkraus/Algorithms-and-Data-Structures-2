/**
 * An interface for a codebook to be used for LZW expansion */

public interface ExpansionCodeBookInterface {

   /**
     * 
     * @return the logical size of the codebook
     */
    public int size();
  
    /**
     * @param flushIfFull a boolean indicating whether the codebook should be 
     *                    flushed when it is at its maximum capacity
     * @return the effective value of the codeword width
     */
    public int getCodewordWidth(boolean flushIfFull);
  
    
    /**
     * Adds a new string to the codebook in O(1) amortized time
     * @param str the String to be added to the codebook
     * @param flushIfFull a boolean indicating whether the codebook should be 
     *                    flushed when it is at its maximum capacity
     */
    public void add(String str, boolean flushIfFull);
  
      
    /**
     * retrieves the string corresponding to a given codeword in O(1) time.
     * @param codeword the int codeword
     * @return the string corresponding to codeword in the codebook; returns null
     * if codeword doesn't exist in the codebook
     */
    public String getString(int codeword);
  
  }