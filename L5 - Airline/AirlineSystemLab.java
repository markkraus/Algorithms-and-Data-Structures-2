/**
 * AirlineSystem.java
 *
 * Implements a Digraph using Adjacency Lists to parse and print graph data
 * containing airline route information.
 *
 * javadoc compilation command:
 * $ javadoc -verbose -private -d ./docs -version -author ./AirlineSystem.java
 *
 * java compilation/run command:
 * $ javac ./AirlineSystem.java
 * $ java AirlineSystem
 *
 * @author  Sherif Khattab
 * @version 1.0
 * @since   2023-03-16
 */

import java.util.*;
import java.io.*;

public class AirlineSystemLab {

	private String [] cityNames;
	private Digraph G;

	/**
	 * Parses the input file into an AirlineSystem object, and output the
	 * Digraph Adjacency List to the console.
	 *
	 * @param args program execution runtime argument string array.
	 */
	public static void main (String[] args) throws IOException {
		AirlineSystemLab airline = new AirlineSystemLab();
		airline.readGraph();
		airline.printGraph();
	}

	/**
	 * Reads in an airline route information data file and parses the contents
	 * to build a representing Digraph. The user input will be a relative
	 * filepath to the data file, which can be opened with a file read stream.
	 */
	public void readGraph () throws IOException {
		Scanner inScan = new Scanner(System.in);
		System.out.println("Please enter graph filename:");
		String fileName = inScan.nextLine();
		Scanner fileScan = new Scanner(new FileInputStream(fileName));

		//TODO1: Complete this method to read the graph from the input file
		//      (e.g., data1.txt)

    int numCities = Integer.parseInt(fileScan.nextLine()); // Get the number of cities in the file

    // Read each city and add them to the array
    cityNames = new String[numCities]; 
    for (int i = 0; i < numCities; i++) {
      cityNames[i] = fileScan.nextLine();
    }

    // Make the new Digraph
    G = new Digraph(numCities);
    while (fileScan.hasNextLine()) {
      // Read connections and add them to the graph
      String line = fileScan.nextLine();
      String[] parts = line.split("\\s+"); // Splitting the source and destination to be added as edge

      // Adjust for zero-indexing
      int source = Integer.parseInt(parts[0]) - 1;
      int destination = Integer.parseInt(parts[1]) - 1;

      // Add the edge to the graph
      G.addEdge(new DirectedEdge(source, destination));
    }

		fileScan.close();
	}

	/**
	 * Prints the Adjacency List of the Digraph G to the console. The neighbor
	 * list will be printed with a whitespace separator, and the neighbor list
	 * will be printed for all vertices, even if a vertex does not have any
	 * adjacent neighbor city.
	 */
	public void printGraph () {

		//TODO2: Complete this method to print the graph as in the expected output
		//      files (e.g., output1.txt)

    for (int v = 0; v < G.v; v++) {
      // Print the source city
      System.out.print(cityNames[v] + ": ");
      for (DirectedEdge e : G.adj(v)) {
        // Print each destination to be a part of the list if it exists
        System.out.print(cityNames[e.to()] + " ");
      }
      // Add a new line in between each city
      System.out.println();
    }
	}

	/**
	*  The <tt>Digraph</tt> class represents an directed graph of vertices
	*  named 0 through v-1. It supports the following operations: add an edge to
	*  the graph, iterate over all of edges leaving a vertex.Self-loops are
	*  permitted.
	*/
	private class Digraph {

		private final int v;
		private int e;
		private LinkedList<DirectedEdge>[] adj;

		/**
		 * Create an empty digraph with v vertices.
		 *
		 * @param v The number of vertices of the Digraph.
		 */
		public Digraph (int v) {
			if (v < 0) throw new RuntimeException("Number of vertices must be nonnegative");
			this.v = v;
			this.e = 0;
			@SuppressWarnings("unchecked")
			LinkedList<DirectedEdge>[] temp =
			(LinkedList<DirectedEdge>[]) new LinkedList[v];
			adj = temp;
			for (int i = 0; i < v; i++) adj[i] = new LinkedList<DirectedEdge>();
		}

		/**
		 * Add the edge e to this digraph.
		 *
		 * @param edge A directed edge between two existing vertices in the Digraph.
		 */
		public void addEdge (DirectedEdge edge) {
			int from = edge.from();
			adj[from].add(edge);
			e++;
		}


		/**
		 * Return the edges leaving vertex v as an Iterable.
		 * To iterate over the edges leaving vertex v, use foreach notation:
		 * <tt>for (DirectedEdge e : graph.adj(v))</tt>.
		 *
		 * @param v Vertex id
		 * @return A DirectedEdge Iterable Object
		 */
		public Iterable<DirectedEdge> adj (int v) {
			return adj[v];
		}
	}


	/**
	 *  The DirectedEdge class represents an edge in an directed graph.
	 */
	private class DirectedEdge {
		private final int u;
		private final int v;

		/**
		 * Create a directed edge from vertex u to v.
		 *
		 * @param u Source vertex id
		 * @param v Destination vertex id
		 */
		public DirectedEdge (int u, int v) {
			this.u = u;
			this.v = v;
		}

		/**
		 * Returns the source vertex of the edge.
		 *
		 * @return vertex id of source vertex.
		 */
		public int from () {
			return u;
		}

		/**
		 * Returns the destination vertex of the edge.
		 *
		 * @return vertex id of destination vertex.
		 */
		public int to () {
			return v;
		}
	}
}
