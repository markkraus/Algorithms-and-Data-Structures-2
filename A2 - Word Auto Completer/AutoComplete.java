/**
 * An implementation of the AutoCompleteInterface using a DLB Trie.
 * @author Mark Kraus
 */


 public class AutoComplete implements AutoCompleteInterface {

  private DLBNode root;                 // Points to the first character added to the trie
  private StringBuilder currentPrefix;  // Stores the user inputted string through the course of the program
  private DLBNode currentNode;          // Pointer to the node being worked on throughout the program

  public AutoComplete(){
    root = null;
    currentPrefix = new StringBuilder();
    currentNode = null;
  }

  /**
   * Adds a word to the dictionary in O(alphabet size*word.length()) time
   * @param word the String to be added to the dictionary
   * @return true if add is successful, false if word already exists
   * @throws IllegalArgumentException if word is the empty string
   */
    public boolean add(String word) {
    if (word.isEmpty()) {
      // String is empty - invalid
      throw new IllegalArgumentException("Word cannot be empty");
    }

    // Get the first letter
    char c = word.charAt(0);

    // Start from the root
    currentNode = root;

    /* These conditions point our currentNode
     to the correct char value of the root of 
     the trie if we don't already point to it */
    if (currentNode == null) {
      // 1. Trie is empty; first addition to the trie
      currentNode = new DLBNode(c);
      root = currentNode;
    } else if (currentNode.data != c) {
      // 2. Root character is different from root
      currentNode = findSiblingNode(currentNode, c);
      if (currentNode == null || currentNode.data != c) {
        // Did not find the first letter to be a root sibling - add it as one
        currentNode = addRootNode(c);
      }
    } 

    /* For each character c in word,
     add it to the trie if it does not exist already */
    for (int i = 1; i < word.length(); i++) {
      c = word.charAt(i); // get character in word
      currentNode = findOrAddChildNode(currentNode, c); // find its node or create it
    }

    if (!currentNode.isWord) {
      // Word DNE in the trie yet - add it
      currentNode.isWord = true; // make current node a word node
      increaseWordCount(currentNode); // Recursively add to each node's word count up the trie
      return true;
    }

    // Current string is already a word
    return false;
  }

  /**
   * Increases each node's word count size recursively in the subtree rooted at a root letter.
   * Works from the leaf nodelet up to the root of the subtree
   * 
   * @param node: The current nodelet being added to
   */
  private void increaseWordCount(DLBNode node) {
    if (node == null) {
      // The root was passed
      return;
    }

    // Add a word count to the current node
    node.size++;
    increaseWordCount(node.parent); // Go the previous letter in the word
  }

  /**
   * Find the sibling of 'node' that contains the data 'c'.
   * @param node: Node where search begins
   * @param c: the desired data of the node we are searching for
   * @return the found node of data 'c', or null if it DNE
   */
  private DLBNode findSiblingNode (DLBNode node, char c) {

    if(node == null) {
      // Passed-in node already DNE, so no siblings exist
      return null;
    }

    while(node != null) {
      // Search all existing siblings
      if(node.data == c) {
        // We found the node being searched for
        return node;
      } else {
        // We haven't found the node yet - go to the next sibling
        node = node.nextSibling;
      }
    }

    // The node doesn't exist in the sibling list
    return null;
  }

  /**
    * Adds a new root sibling to the DLB Trie.
    * @param firstLetter: root character of a word
    * @return the newly added root sibling nodelet
    */
    private DLBNode addRootNode(char firstLetter) {
      DLBNode newRoot = root; // Point to the root of the trie
      DLBNode prevNode;

      while(newRoot.nextSibling != null) {
        // Go to the end of the root sibling list
        newRoot = newRoot.nextSibling;
      }

      // Add the new node as the last sibling in the sibling list
      prevNode = newRoot;
      newRoot = new DLBNode(firstLetter);
      prevNode.nextSibling = newRoot;
      newRoot.previousSibling = prevNode;
      
      return newRoot;
   }

   /**
    *  Adds a child of 'parent', and updates sibling relationships accordingly,
    *  or returns an already existing child of 'data'.
    * @param parent: The parent node to append the child to
    * @param data: The data of the child node to find/add
    * @return the child node of 'data' with updated child properties/relations
    */
  private DLBNode findOrAddChildNode(DLBNode parent, char data) {
    DLBNode child = parent.child;   // Link to the child of the passed-in node
    DLBNode previousChild = null; // Used for creating a connection amongst sibling nodes

    while(child != null && child.data != data) {
      // Check if the passed-in data already exists as a child of the parent node
      previousChild = child;
      child = child.nextSibling;
    }

    if (child == null) {
      // No child of 'data' existed - create a new child of 'data'
      child = new DLBNode(data);

      if(previousChild == null) {
        // There was no preexisting child node before - make 'data' the first child
        parent.child = child;
      } else {
        // There was already a child - make the newly created child a sibling of the previous child
        previousChild.nextSibling = child;
        child.previousSibling = previousChild;
      }

      // Make the newly created child's parent node the passed-in node
      child.parent = parent;
    } 
    
    return child;
  }


  /**
   * appends the character c to the current prefix in O(alphabet size) time. 
   * This method doesn't modify the dictionary.
   * @param c: the character to append
   * @return true if the current prefix after appending c is a prefix to a word 
   * in the dictionary and false otherwise
   */
    public boolean advance(char c){
      currentPrefix.append(c); // Add the user input to the prefix

      // Find the node path of the prefix
      currentNode = getNode(root, currentPrefix.toString(), 0);

      if(currentNode == null) {
        // Path DNE, so prefix is not valid
        return false;
      }

      // Every node of the prefix string was found in the trie
      return true;
    }


  /**
   * removes the last character from the current prefix in O(alphabet size) time. This 
   * method doesn't modify the dictionary.
   * @throws IllegalStateException if the current prefix is the empty string
   */
    public void retreat() {
    if (currentPrefix.length() == 0) {
      // Tried to delete from the prefix when it was already empty
      throw new IllegalStateException();
    }

    // Delete the last letter in the prefix
    currentPrefix.deleteCharAt(currentPrefix.length() - 1);

    if (currentPrefix.length() > 0) {
      // Still have letters in the prefix, so find the last letter of the prefix
      currentNode = getNode(root, currentPrefix.toString(), 0);
    }
  }

  /**
   * resets the current prefix to the empty string in O(1) time
   */
    public void reset(){
      // Make a brand new StringBuilder that wipes the old data
      currentPrefix = new StringBuilder();
    }
    
  /**
   * @return true if the current prefix is a word in the dictionary and false
   * otherwise. The running time is O(1).
   */
    public boolean isWord(){
      // We don't have a node, or it DNE as a word
      if(currentNode == null || !currentNode.isWord) return false;

      // Node is a word in the trie
      return true;
    }

  /**
   * adds the current prefix as a word to the dictionary (if not already a word)
   * The running time is O(alphabet size*length of the current prefix). 
   */
    public void add(){
      if (currentNode == null || !currentNode.isWord) {
        // Node DNE, or is not a word in the trie - Make the current sequence a word in the trie
        String newWord = currentPrefix.toString();
        add(newWord);
      }
    }

  /** 
   * @return the number of words in the dictionary that start with the current 
   * prefix (including the current prefix if it is a word). The running time is 
   * O(1).
   */
    public int getNumberOfPredictions(){
      // Node DNE
      if (currentNode == null) return 0;

      // Node exists and has a word count
      return currentNode.size;
    }
  
  /**
   * retrieves one word prediction for the current prefix. The running time is 
   * O(prediction.length())
   * @return a String or null if no predictions exist for the current prefix
   */
    public String retrievePrediction() {
      // Holds the word we will predict based off of the current prefix
      String prediction = currentPrefix.toString();

      while (currentNode != null) {
        if (currentNode.isWord) {
          // We have a prediction with the prefix - get the last node of that word
          currentNode = getNode(root, currentPrefix.toString(), 0);

          // Send that prediction back
          return prediction;
        } 
        // We do not yet have a prediction from the prefix - go down the trie, adding the letter along the way
        currentNode = currentNode.child;
        prediction += currentNode.data;
      }
    
      // We couldn't find a prediction
      return null;
    }

  /* ==============================
   * Helper methods for debugging.
   * ==============================
   */

  //print the subtrie rooted at the node at the end of the start String
  public void printTrie(String start){
    System.out.println("==================== START: DLB Trie Starting from \""+ start + "\" ====================");
    if(start.equals("")){
      printTrie(root, 0);
    } else {
      DLBNode startNode = getNode(root, start, 0);
      if(startNode != null){
        printTrie(startNode.child, 0);
      }
    }
    
    System.out.println("==================== END: DLB Trie Starting from \""+ start + "\" ====================");
  }

  //a helper method for printTrie
  private void printTrie(DLBNode node, int depth){
    if(node != null){
      for(int i=0; i<depth; i++){
        System.out.print(" ");
      }
      System.out.print(node.data);
      if(node.isWord){
        System.out.print(" *");
      }
      System.out.println(" (" + node.size + ")");
      printTrie(node.child, depth+1);
      printTrie(node.nextSibling, depth);
    }
  }

  //return a pointer to the node at the end of the start String 
  //in O(start.length() - index)
  private DLBNode getNode(DLBNode node, String start, int index){
    if(start.length() == 0){
      return node;
    }
    DLBNode result = node;
    if(node != null){
      if((index < start.length()-1) && (node.data == start.charAt(index))) {
          result = getNode(node.child, start, index+1);
      } else if((index == start.length()-1) && (node.data == start.charAt(index))) {
          result = node;
      } else {
          result = getNode(node.nextSibling, start, index);
      }
    }
    return result;
  } 

  //The DLB node class
  private class DLBNode{
    private char data;
    private int size;
    private boolean isWord;
    private DLBNode nextSibling;
    @SuppressWarnings("unused")
    private DLBNode previousSibling;
    private DLBNode child;
    private DLBNode parent;

    private DLBNode(char data){
        this.data = data;
        size = 0;
        isWord = false;
        nextSibling = previousSibling = child = parent = null;
    }
  }
}