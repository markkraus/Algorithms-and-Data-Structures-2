import java.util.*;
import java.io.*;

public class AirlineSystem implements AirlineInterface {
  private String[] cityNames; // Holds the name of each city in the input file
  private Digraph G;          // Graph object containing city vertices as well as edge's source city, destination city, price, and distance

  public boolean loadRoutes(String fileName) {
    try {
      // Read the file
      BufferedReader bf = new BufferedReader(new FileReader(fileName));

      // Map each city to a designated index
      int numCities = Integer.parseInt(bf.readLine()); // Get the total amount of cities in the graph
      cityNames = new String[numCities];
      for (int i = 0; i < numCities; i++) {
        // Read each city from the input file and add it
        String city = bf.readLine();
        cityNames[i] = city;
      }

      // Create a new Digraph with vertices for each city
      G = new Digraph(numCities);

      // Read and store flight values
      while (bf.ready()) {
        String line = bf.readLine();
        String[] flightInfo = line.split("\\s+");

        // Gather each piece of data from the flight
        int source = Integer.parseInt(flightInfo[0]) - 1; // Subtract for zero-indexing
        int destination = Integer.parseInt(flightInfo[1]) - 1;
        int distance = Integer.parseInt(flightInfo[2]);
        double price = Double.parseDouble(flightInfo[3]);

        // Add each flight between the source/destination vertices
        G.addEdge(new DirectedEdge(source, destination, distance, price));
        G.addEdge(new DirectedEdge(destination, source, distance, price));
      }

      // No issues occured
      bf.close();
      return true;
    }
    catch (IOException e) {
      // Some error occured
      e.printStackTrace();
    }
    return false;
  }

  public Set<String> retrieveCityNames() {
    Set<String> cities = new HashSet<String>();

    // Add every city we read from the file to the set
    cities.addAll(Arrays.asList(cityNames));

    return cities;
  }

  public Set<Route> retrieveDirectRoutesFrom(String city) throws CityNotFoundException {
    Set<Route> flights = new HashSet<>();

    // Find the passed-in city, or throw exception if it's unknown
    int cityIndex = findCityIndex(city);
    if (cityIndex == -1) throw new CityNotFoundException(city);
    
    // Add each flight from the passed-in city
    for (DirectedEdge edge : G.adj[cityIndex]) {
      flights.add(new Route(cityNames[edge.from()], cityNames[edge.to()], edge.distance(), edge.price()));
    }

    return flights;
  }

  public Set<ArrayList<String>> fewestStopsItinerary(String source, String destination) throws CityNotFoundException {
    Set<ArrayList<String>> fewestStopsPath = new HashSet<>();

    // Find the passed-in cities, or throw exception if they're unknown
    int sourceIndex = findCityIndex(source);
    int destinationIndex = findCityIndex(destination);
    if (sourceIndex == -1 || destinationIndex == -1) throw new CityNotFoundException(source + " " + destination);

    // Find the shortest path using Breadth-First Search
    G.bfs(sourceIndex);

    // Find every shortest path from destination - source
    for (DirectedEdge edge : G.adj[destinationIndex]) {
      // Check each edge coming from the destination vertex
      if (G.hops[destinationIndex] - G.hops[edge.to()] == 1) {
        // Neighbor vertex will be 1 less hop from the source - it's a shortest route
        ArrayList<String> path = new ArrayList<>();

        // Add the destination vertex to the path and move to the neighbor
        path.add(0, cityNames[destinationIndex]);
        int v = edge.to();
        while (v != sourceIndex) {
          // Add each vertex along the path to the source
          path.add(0, cityNames[v]);
          v = G.edgeTo[v];
        }
        // Add the source vertex to the path and add the path to the Set
        path.add(0, source); 
        fewestStopsPath.add(path);
      }
    }

    return fewestStopsPath;
  }

  public Set<ArrayList<Route>> shortestDistanceItinerary(String source, String destination) throws CityNotFoundException {
    Set<ArrayList<Route>> shortestDistancePaths = new HashSet<>();

    // Find the passed-in cities, or throw exception if they're unknown
    int sourceIndex = findCityIndex(source);
    int destinationIndex = findCityIndex(destination);
    if (sourceIndex == -1 || destinationIndex == -1) throw new CityNotFoundException(source + " " + destination);

    // Find the shortest path based on distance
    G.dijkstras(sourceIndex, destinationIndex);

    // Find all shortest paths with the same minimum distance
    int shortestDistance = G.distance[destinationIndex];
    Deque<Integer> stack = new ArrayDeque<>();
    stack.push(destinationIndex);
    findAllShortestPaths(sourceIndex, destinationIndex, shortestDistance, stack, shortestDistancePaths);

    return shortestDistancePaths;
  }

  public Set<ArrayList<Route>> cheapestItinerary(String source, String destination) throws CityNotFoundException {
    Set<ArrayList<Route>> cheapestPaths = new HashSet<>();

    // Find the passed-in cities, or throw an exception if they're unknown
    int sourceIndex = findCityIndex(source);
    int destinationIndex = findCityIndex(destination);
    if (sourceIndex == -1 || destinationIndex == -1) throw new CityNotFoundException(source + " " + destination);

    // Find the cheapest path
    G.dijkstrasPrice(sourceIndex, destinationIndex);


    // Traverse all edges and find paths with the same minimum cost
    double cheapestPrice = G.price[destinationIndex];
    Deque<Integer> stack = new ArrayDeque<>();
    G.visited = new boolean[G.v];
    stack.push(destinationIndex);
    findAllCheapestPaths(sourceIndex, destinationIndex, cheapestPrice, stack, cheapestPaths, G.visited);

    return cheapestPaths;
  }

  public Set<ArrayList<Route>> cheapestItinerary(String source, String transit, String destination) throws CityNotFoundException {
    Set<ArrayList<Route>> cheapestPaths = new HashSet<>();

    // Make sure all cities exist, or throw an error if any of them are unknown
    int sourceIndex = findCityIndex(source);
    int transitIndex = findCityIndex(transit);
    int destinationIndex = findCityIndex(destination); 
    if (sourceIndex == -1 || transitIndex == -1 || destinationIndex == -1) throw new CityNotFoundException(source + " " + transit + " " + destination);

    // Get price of source -> transit
    G.dijkstrasPrice(sourceIndex, transitIndex);
    double toTransitPrice = G.price[transitIndex];

    // Get price of transit -> destination
    G.dijkstrasPrice(transitIndex, destinationIndex);
    double fromTransitPrice = G.price[destinationIndex];

    // Sets to hold all possible paths from source -> transit -> destination with cheapest price
    Set<ArrayList<Route>> sourceToTransitPaths = new HashSet<>();
    Set<ArrayList<Route>> transitToDestinationPaths = new HashSet<>();

    // Find all cheapest paths from source -> transit
    Deque<Integer> stack = new ArrayDeque<>();
    G.visited = new boolean[G.v];
    stack.push(transitIndex);
    findAllCheapestPaths(sourceIndex, transitIndex, toTransitPrice, stack, sourceToTransitPaths, G.visited);

    stack.clear();

    // Find all cheapest paths from transit -> destination
    G.visited = new boolean[G.v];
    stack.push(destinationIndex);
    findAllCheapestPaths(transitIndex, destinationIndex, fromTransitPrice, stack, transitToDestinationPaths, G.visited);

    // Combine the paths to form the complete itineraries
    for (ArrayList<Route> sourceToTransitPath : sourceToTransitPaths) {
      for (ArrayList<Route> transitToDestinationPath : transitToDestinationPaths) {
        ArrayList<Route> completePath = new ArrayList<>(sourceToTransitPath);
        completePath.addAll(transitToDestinationPath.subList(0, transitToDestinationPath.size()));
        cheapestPaths.add(completePath);
      }
    }
    return cheapestPaths;
  }

  public Set<Set<Route>> getMSTs() {
    Set<Set<Route>> minimumSpanningTrees = new HashSet<>();

    // Keep track of visited vertices
    boolean[] visited = new boolean[G.v];

    // Find the MST that accounts for every vertex
    for (int i = 0; i < G.v; i++) {
      if (!visited[i]) {
        Set<Route> minimumSpanningTree = new HashSet<>();
        Set<DirectedEdge> mstEdges = new HashSet<>();
        PriorityQueue<DirectedEdge> pq = new PriorityQueue<>(
            Comparator.comparingInt(DirectedEdge::distance));
        visited[i] = true;

        // Add all edges of the current vertex to the priority queue
        for (DirectedEdge edge : G.adj(i)) {
          pq.add(edge);
        }

        // Prim's algorithm to find the minimum spanning tree
        while (!pq.isEmpty()) {
          DirectedEdge currentEdge = pq.poll();
          int u = currentEdge.from();
          int v = currentEdge.to();
          if (visited[u] && visited[v]) {
            continue;
          }

          mstEdges.add(currentEdge); // Add the current edge to the MST

          // Add neighboring edges to the priority queue and mark vertices as visited
          if (!visited[u]) {
            visited[u] = true;
            for (DirectedEdge edge : G.adj(u)) {
              pq.add(edge);
            }
          }
          if (!visited[v]) {
            visited[v] = true;
            for (DirectedEdge edge : G.adj(v)) {
              pq.add(edge);
            }
          }
        }

        // Convert edges to routes and add the MST to the set of MSTs
        for (DirectedEdge edge : mstEdges) {
          minimumSpanningTree.add(new Route(cityNames[edge.from()], cityNames[edge.to()], edge.distance(), edge.price()));
        }
        minimumSpanningTrees.add(minimumSpanningTree);
      }
    }
    return minimumSpanningTrees;
  }

  public Set<ArrayList<Route>> tripsWithin(String city, double budget) throws CityNotFoundException {
    Set<ArrayList<Route>> validPaths = new HashSet<>();

    // Find the city index, or throw an exception if it's unknown
    int sourceIndex = findCityIndex(city);
    if (sourceIndex == -1) throw new CityNotFoundException(city);

    // Initialize the current path price, and set the source city as visited to start
    double pathPrice = 0;
    boolean[] visited = new boolean[G.v];
    visited[sourceIndex] = true;

    // Create an initial path with the starting city
    ArrayList<Route> initialPath = new ArrayList<>();

    // Find all paths from souce city that are under the budget
    tripsWithinRecursive(sourceIndex, budget, pathPrice, visited, initialPath, validPaths);

    return validPaths;
  }

  public Set<ArrayList<Route>> tripsWithin(double budget) {
    Set<ArrayList<Route>> validPaths = new HashSet<>();

    // Initialize the current path price
    double pathPrice = 0;

    // Create an initial path
    ArrayList<Route> initialPath = new ArrayList<>();

    // Iterate over all vertices to start the recursive search from each city
    for (int sourceIndex = 0; sourceIndex < G.v; sourceIndex++) {
      // Reset visited array for each source city
      boolean[] visited = new boolean[G.v];
      visited[sourceIndex] = true;

      // Find all paths from the source city that are under the budget
      tripsWithinRecursive(sourceIndex, budget, pathPrice, visited, initialPath, validPaths);
    }

    return validPaths;
  }

  /**
   * Finds the index of city in the cityNames array.
   * @param city the desired city to find.
   * @return the index of the city, or -1 if it's not found.
   */
  private int findCityIndex(String city) {
    int cityIndex = -1;

    // Find the passed-in city, or throw an exception if it's not found
    for (int i = 0; i < cityNames.length; i++) {
      if (city.equalsIgnoreCase(cityNames[i])) {
        // City was found
        cityIndex = i;
        break;
      }
    }

    return cityIndex;
  }

  /**
   * Recursively finds all the shortest paths between the source and the
   * destination based on distance.
   *
   * @param sourceIndex      Index of the source city.
   * @param destinationIndex Index of the destination city.
   * @param shortestDistance Shortest distance between the source and the current
   *                         destination.
   * @param stack            Stack used for tracking the current path.
   * @param paths            Set containing all the paths under the shortest
   *                         distance.
   */
  private void findAllShortestPaths(int sourceIndex, int destinationIndex, int shortestDistance, Deque<Integer> stack,
      Set<ArrayList<Route>> paths) {
    // If the shortest distance is 0, the destination has been reached, and the
    // complete path is added to the set
    if (shortestDistance == 0) {
      // Destination has been reached with the shortest distance; add the path
      ArrayList<Route> path = new ArrayList<>();
      int prev = -1;
      // Iterate through the stack to create Route objects and add them to the path
      for (int v : stack) {
        if (prev != -1) {
          path.add(new Route(cityNames[prev], cityNames[v], G.distance[v] - G.distance[prev], G.price[v] - G.price[prev]));
        }
        prev = v;
      }
      paths.add(path);
    } else {
      // Traverse all the edges connected to the current destination vertex
      for (DirectedEdge edge : G.adj[destinationIndex]) {
        // Make sure this edge will be short enough for now
        if (G.distance[destinationIndex] - G.distance[edge.to()] == edge.distance()) {
          // Push the edge's destination to the stack and recursively find all shortest
          // paths from that point
          stack.push(edge.to());
          findAllShortestPaths(sourceIndex, edge.to(), shortestDistance - edge.distance(), stack, paths);
          // Remove the last element from the stack to backtrack and explore other paths
          stack.pop();
        }
      }
    }
  }

  /**
   * Recursively finds all the cheapest paths between the source and the 
   * destinaton based on distance.
   * 
   * @param sourceIndex      Index of the source city.
   * @param destinationIndex Index of the destination city.
   * @param cheapestPrice    Cheapest price between the source and the current
   *                         destination.
   * @param stack            Stack used for tracking the current path.
   * @param paths            Set containing all the paths under the cheapest
   *                         price.
   * @param visited          Tracks what vertices are in the current path.
   * 
   */
  private void findAllCheapestPaths(int sourceIndex, int destinationIndex, double cheapestPrice, Deque<Integer> stack, Set<ArrayList<Route>> paths, boolean[] visited) {
    if (cheapestPrice == 0 && destinationIndex == sourceIndex) {
      // We've hit the source vertex, and our total is equal to the cheapest itinerary - add the path
      ArrayList<Route> path = new ArrayList<>();
      int prev = -1;

      // Iterate through the stack to create Route objects and add them to the path
      for (int v : stack) {
        if (prev != -1) {
          path.add(new Route(cityNames[prev], cityNames[v], Math.abs(G.distance[v] - G.distance[prev]), Math.abs(G.price[v] - G.price[prev])));
        }
        prev = v;
      }

      // Add the completed path to the set of paths
      paths.add(path);
      visited[destinationIndex] = false; // Set false for recursive backtracking
    } else {
      // Haven't finished the path yet - check edge connecting edge to the vertex to see if
      // we can possibly go down that path
      for (DirectedEdge edge : G.adj[destinationIndex]) {
        visited[destinationIndex] = true; // mark the current vertex as visited
        int nextVertex = edge.to(); // get a connecting edge
        if (!visited[edge.to()] && (Math.abs((int) G.price[destinationIndex] - (int) G.price[nextVertex]) <= cheapestPrice)) {
          // This path is cheap enough to go down for now

          // Push the connecting vertex to the stack and recursively find all cheapest paths
          // from that point
          visited[nextVertex] = true;
          stack.push(nextVertex);
          double originalCheapestPrice = cheapestPrice;
          cheapestPrice -= edge.price();
          findAllCheapestPaths(sourceIndex, nextVertex, cheapestPrice, stack, paths, visited);

          // We've returned from an unsuccessful path, pop the current vertex off
          stack.pop();

          // Restore the original cheapestPrice after coming back from the recursive call
          cheapestPrice = originalCheapestPrice;

        }
      }
      // Backtracking; mark as unvisited since we are trying a new path
      visited[destinationIndex] = false;
    }
  }

  /**
   * Recursively finds each route less than or equal to a given budget.
   * 
   * @param currentVertex Source city for each call.
   * @param budget        User requested budget.
   * @param pathPrice     The current path price.
   * @param visited       True for vertices already added to the path;
   *                      false otherwise.
   * @param currentPath   The current path of routes.
   * @param validPaths    Set of all paths less than or equal to budget.
   * 
   */
  private void tripsWithinRecursive(int currentVertex, double budget, double pathPrice, boolean[] visited, ArrayList<Route> currentPath, Set<ArrayList<Route>> validPaths) {
    // Check all flights from the current vertex
    for (DirectedEdge edge : G.adj[currentVertex]) {
      // If the city hasn't been visited and the cost is within the budget
      if (!visited[edge.to()] && pathPrice + edge.price() <= budget) {
        // City hasn't been visited and this edge is under the budget

        // mark the city as visited, and update current path price
        visited[edge.to()] = true;
        pathPrice += edge.price();

        // Add the flight to the current path
        currentPath.add(new Route(cityNames[currentVertex], cityNames[edge.to()], edge.distance(), edge.price()));

        // Recursively call to explore additional flights
        tripsWithinRecursive(edge.to(), budget, pathPrice, visited, currentPath, validPaths);

        // Backtrack: Remove the flight and update the path, city visitation, and cost
        visited[edge.to()] = false;
        pathPrice -= edge.price();
        currentPath.remove(currentPath.size() - 1);
      }
    }

    // Add the current path to the set of valid paths
    if (!currentPath.isEmpty()) {
      validPaths.add(new ArrayList<>(currentPath));
    } 
  }

  /**
   * The <tt>Digraph</tt> class represents an directed graph of vertices
   * named 0 through v-1. It supports the following operations: add an edge to
   * the graph, iterate over all of edges leaving a vertex.Self-loops are
   * permitted.
   */
  private class Digraph {

    private final int v;
    private int e;
    private LinkedList<DirectedEdge>[] adj;
    int[] hops;
    private int[] distance; // Number of hops from source vertex to destination vertex
    private double[] price;
    private int[] edgeTo; // Holds the vertex index of the vertex coming to
    private boolean[] visited;

    /**
     * Create an empty digraph with v vertices.
     *
     * @param v The number of vertices of the Digraph.
     */
    private Digraph(int v) {
      if (v < 0)
        throw new RuntimeException("Number of vertices must be nonnegative");
      this.v = v;
      this.e = 0;
      @SuppressWarnings("unchecked")
      LinkedList<DirectedEdge>[] temp = (LinkedList<DirectedEdge>[]) new LinkedList[v];
      adj = temp;
      for (int i = 0; i < v; i++)
        adj[i] = new LinkedList<DirectedEdge>();
    }

    /**
     * Add the edge e to this digraph.
     *
     * @param edge A directed edge between two existing vertices in the Digraph.
     */
    private void addEdge(DirectedEdge edge) {
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
    private Iterable<DirectedEdge> adj(int v) {
      return adj[v];
    }

    private void bfs(int source) {
      Queue<Integer> q = new LinkedList<>();
      boolean[] visited = new boolean[v];
      boolean[] seen = new boolean[v];
      int[] parent = new int[v];
      edgeTo = new int[v];
      hops = new int[v];

      // Add the source city to start the search
      q.add(source);
      seen[source] = true;
      while (!q.isEmpty()) {
        // Remove the head, and mark it as visited
        int w = q.poll();
        visited[w] = true;

        // Add each unseen neighbor of the just-removed vertex
        for (DirectedEdge edge : G.adj[w]) {
          if (!seen[edge.to()]) {
            seen[edge.to()] = true;
            parent[edge.to()] = w;
            hops[edge.to()] = hops[w] + 1; // It took another hop to get to this index
            edgeTo[edge.to()] = w; // The just-deleted node connected to this index
            q.add(edge.to());
          }
        }
      }
    }

    private void dijkstras(int source, int destination) {
      // Initialize data structures for tracking distances and visited vertices
      boolean[] visited = new boolean[v]; // to keep track of visited vertices
      distance = new int[v]; // to keep track of the distance from the source
      edgeTo = new int[v];
      price = new double[v];
      int INFINITY = Integer.MAX_VALUE;

      // Set initial distances to infinity and source distance to 0
      Arrays.fill(distance, INFINITY);
      distance[source] = 0;

      // Create a priority queue to store vertices with the minimum distance
      PriorityQueue<Integer> pq = new PriorityQueue<>(this.v, Comparator.comparingInt(vertex -> distance[vertex]));

      // Add the source vertex to the priority queue
      pq.add(source);

      while (!pq.isEmpty()) {
        // Extract the vertex with the minimum distance
        int u = pq.poll(); // u = source
        visited[u] = true;

        // Traverse the neighboring vertices of the current vertex
        for (DirectedEdge edge : adj(u)) {
          int v = edge.to(); // v = destination
          int weight = edge.distance();
          double priceWeight = edge.price();

          // If the distance to the neighboring vertex can be reduced, update the distance
          // and add it to the priority queue
          if (!visited[v] && distance[u] + weight < distance[v]) {
            distance[v] = distance[u] + weight;
            price[v] = price[u] + priceWeight;
            edgeTo[v] = u;
            pq.add(v);
          }
        }
      }
    }

    private void dijkstrasPrice(int source, int destination) {
      boolean[] visited = new boolean[v]; // to keep track of visited vertices
      distance = new int[v];
      edgeTo = new int[v];
      price = new double[v];
      int INFINITY = Integer.MAX_VALUE;

      // Set initial prices to infinity and source distance to 0
      Arrays.fill(price, INFINITY);
      price[source] = 0;

      // Create a priority queue to store vertices with the minimum distance
      PriorityQueue<Integer> pq = new PriorityQueue<>(this.v, Comparator.comparingInt(vertex -> (int) price[vertex]));

      // Add the source vertex to the priority queue
      pq.add(source);

      while (!pq.isEmpty()) {
        // Extract the vertex with the minimum distance
        int u = pq.poll(); // u = source
        visited[u] = true;

        // Traverse the neighboring vertices of the current vertex
        for (DirectedEdge edge : adj(u)) {
          int v = edge.to(); // v = destination
          double weight = edge.price();

          // If the cost to the neighboring vertex can be reduced, update the distance
          // and add it to the priority queue
          if (!visited[v] && price[u] + weight < price[v]) {
            price[v] = price[u] + weight;
            distance[v] = distance[u] + edge.distance();
            edgeTo[v] = u;
            pq.add(v);
          }
        }
      }
    }
  }

  /**
   * The DirectedEdge class represents an edge in an directed graph.
   */
  private class DirectedEdge {
    private final int u;
    private final int v;
    private final int distance;
    private final double price;

    /**
     * Create a directed edge from vertex u to v.
     *
     * @param u Source vertex id
     * @param v Destination vertex id
     */
    private DirectedEdge(int u, int v, int distance, double price) {
      this.u = u;
      this.v = v;
      this.distance = distance;
      this.price = price;
    }

    /**
     * Returns the source vertex of the edge.
     *
     * @return vertex id of source vertex.
     */
    private int from() {
      return u;
    }

    /**
     * Returns the destination vertex of the edge.
     *
     * @return vertex id of destination vertex.
     */
    private int to() {
      return v;
    }

    /**
     * Returns the distance of the edge.
     * 
     * @return distance between source/destination vertices.
     */
    private int distance() {
      return distance;
    }

    /**
     * Returns the price of the edge.
     * 
     * @return price between source/destination vertices.
     */
    private double price() {
      return price;
    }
  }
}