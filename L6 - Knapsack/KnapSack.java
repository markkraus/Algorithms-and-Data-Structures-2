import java.io.File;
import java.util.Scanner;
import java.util.ArrayList;
import java.util.Arrays;

public final class KnapSack {

	private Item[] items;
	private int nItems;
	private int capacity;
	private int[][] KS; //the dynamic programming table
	private int[][] decision;//remember the decision that we make at each subproblem
	private ArrayList<Item> selectedItems;

	public static void main(String[] args){
		if(args.length != 1){
			System.out.println("Usage: java KnapSack <filename>");
		} else {
			new KnapSack(args[0]);
		}
	}

	public KnapSack(String filename){
    Scanner scan = null;
		try {
				scan = new Scanner(new File(filename));
				nItems = Integer.parseInt(scan.nextLine());
				items = new Item[nItems];
				for(int i=0; i<nItems; i++){
          items[i] = new Item(scan.nextInt(), scan.nextInt());
          scan.nextLine();
				}
				capacity = Integer.parseInt(scan.nextLine());
				KS = new int[nItems+1][capacity+1];
        decision = new int[nItems+1][capacity+1];
				selectedItems = new ArrayList<>();

		} catch (Exception e) {
			System.out.println("Error reading file " + filename + ": " + e.getMessage());
		} finally{
      scan.close();
    }

		solve();

		System.out.println("The maximum value that can fit with a knapsack of size "
											+ capacity + " is " + KS[nItems][capacity]);

		buildSolution();
		System.out.println("The selected items are: " + selectedItems.toString());

	}

	private void solve(){
    for(int i=0; i<=capacity; i++){
      KS[0][i] = 0;
    }
    for(int i=0; i<=nItems; i++){
      KS[i][0] = 0;
    }
		for(int i=1; i<= nItems; i++){
      for(int j=1; j<= capacity; j++){
          if(j >= items[i-1].weight){
            if(items[i-1].value + KS[i-1][j-items[i-1].weight] > KS[i-1][j]){
              KS[i][j] = items[i-1].value + KS[i-1][j-items[i-1].weight];
              decision[i][j] = 1; //item i-1 taken
            } else {
              //Update KS[i][j] and decision[i][j] to reflect that
              //item i-1 wasn't taken
              KS[i][j] = KS[i-1][j];
              decision[i][j] = 0;
            }
          } else {
						KS[i][j] = KS[i-1][j];
						decision[i][j] = 0; //item i-1 not taken
					}
      }
    }
  }

	private void buildSolution(){
		int curCol = capacity;
    for(int curRow = nItems; curRow >= 0; curRow--){
      if(decision[curRow][curCol] == 1){ //item curRow-1 was taken
        selectedItems.add(items[curRow-1]);
				//Update the current column (curCol) in the KS array
				//Hint: It should move to the left
        curCol -= items[curRow-1].weight;
      }
    }
	}

  private class Item {
    private int weight;
    private int value;
    private Item(int w, int v){
      weight = w;
      value = v;
    }
    public String toString(){
      return "(weight=" + weight + ", value=" + value + ")";
    }
  }

}
