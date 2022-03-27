package main;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import model.DFSCode;
import model.Edge;
import model.Graph;
import model.History;
import model.PDFS;
import model.Projected;
import model.SimpleEdge;
import model.Vertex;

import java.util.NavigableMap;
import java.util.Set;
import java.util.TreeMap;
import java.util.Vector;
import java.util.function.BiFunction;

public class TEDSLightProcessor {
    private ArrayList<Graph> TRANS;
    private DFSCode DFS_CODE;
    private DFSCode DFS_CODE_IS_MIN;
    private Graph GRAPH_IS_MIN;

    private long ID;
    //private long minSup;
    //private long arg.minNodeNum;
    //private long arg.maxNodeNum;
    private boolean directed;
    private FileWriter os;
    private ArrayList<Graph> allGraphs;
    private HashMap<Integer, Set<Integer>> CoveredEdges_OriginalGraphs;
    private Set<Integer>            allCoveredEdges;
    private ArrayList<Set<Integer>> CoveredEdges_EachPattern;
    //// |priv(pattern)|: private edges for each pattern
    private ArrayList<Integer> Priv_pattern; 
    ////            rcov:  reverse coverage for each edge
    private HashMap<Integer, Set<Integer>> Rcov_edge;
    ////            |cov|: number of edges covered
    private Integer numberofcovered;
    ////            C_min: minimum pattern
    private Integer minimumpattern_score;
    private Integer minimumpattern_id;
    ////            rpriv:reverse private edges
    private HashMap<Integer, Set<Integer>> Rpriv_i;
    // Single vertex handling stuff [graph][vertexLabel] = count.
    private NavigableMap<Integer, NavigableMap<Integer, Integer>> singleVertex;
    private NavigableMap<Integer, Integer> singleVertexLabel;
    private Arguments arg;
    
    //////////////// Simple Index////////////////
    ///// [patternid] = set<Integer>, it indicates the edges which are contained by the patternid
    private HashMap<Long, Set<Integer>> CoveredEdges_patterns;
    ///// [edgeid] = set<Integer>, it indicates the patterns that contains the edge with  edgeid
    private HashMap<Integer, Set<Integer>> PatternsID_edges;
    
    public TEDSLightProcessor() {
        TRANS = new ArrayList<>();
        DFS_CODE = new DFSCode();
        DFS_CODE_IS_MIN = new DFSCode();
        GRAPH_IS_MIN = new Graph();
        singleVertex = new TreeMap<>();
        singleVertexLabel = new TreeMap<>();
        /////////////////////////////////////////////////
        //allDFSCodes = new ArrayList<DFSCode>();
        allGraphs = new ArrayList<Graph>();
        allCoveredEdges = new HashSet<Integer>();
        CoveredEdges_OriginalGraphs = new HashMap<Integer, Set<Integer>>();
        Priv_pattern = new ArrayList<Integer>(); 
        Rcov_edge    = new HashMap<Integer, Set<Integer>>();
        numberofcovered = 0;
        minimumpattern_score  = -1;
        Rpriv_i = new HashMap<Integer, Set<Integer>>();
        CoveredEdges_EachPattern = new ArrayList<Set<Integer>>(); 
        
        CoveredEdges_patterns = new HashMap<Long, Set<Integer>>();
        PatternsID_edges= new HashMap<Integer, Set<Integer>>();
    }

    void run(FileReader reader, FileWriter writers, Arguments arguments) throws IOException {
        os = writers;
        ID = 0;
        directed = false;
        arg = arguments;
        if(arg.isLightVersion) {
        	readLight(reader);
        	System.out.println("Number of sampled graphs: " + this.TRANS.size());
        	//return ;
        }else {
        	read(reader);
        }
        
        Long Time1 = System.currentTimeMillis();	
        if(arg.hasInitialPatternGenerator && !arg.strategy.equals("greedy")) {
        	InitialPatternGenerator();
        }
        Long Time2 = System.currentTimeMillis();
		System.out.println("InitialPatternGenerator Time(s) : " + (Time2 - Time1)*1.0 /1000);
		
        runIntern();
        
       // reportIndexSize();
        
        Long Time3 = System.currentTimeMillis();
		System.out.println("Swapping Time(s) : " + (Time3 - Time2)*1.0 /1000);
        
        if(!arg.strategy.equals("greedy")) {
        	  int count = 0;
              for(Graph g : allGraphs) {
              	newreport(g, count++);
              }
              System.out.println("TopK, After Swapping, Number of covered edges: " + allCoveredEdges.size());
              int totalegdes = 0;
              for(int i=0;i< TRANS.size();i++) {
              	totalegdes += TRANS.get(i).getEdgeSize();
               }
              System.out.println("totalegdes : " +  totalegdes);
              System.out.println("Coverage rate : " + allCoveredEdges.size()*1.0 / totalegdes);
        }else {
        	//int count = 0;
        	//for(Graph g : allGraphs) {
            //  	newreport(g, count++);
            //  }
        	
        	Set<Integer> coverededges_curall = new HashSet<Integer>();
        	Set<Long> selectPatternIndex  	 = new HashSet<Long>();
        	int tempcount = 0;
        	System.out.println("allGraphs.size():" + allGraphs.size());
        	while(tempcount < arg.numberofpatterns) {
        		long maxid = -1;
        		int maxgain = -1;
        		for(int k=0;k<allGraphs.size() ;k++) {
        			long id = (long)k;
        			if(selectPatternIndex.contains(id))  continue; 
        			int gain = 0;
        			Set<Integer> converages = CoveredEdges_patterns.get(id);
        			for(Integer a: converages) {
        				if(!coverededges_curall.contains(a)) gain++;
        			}
        			if(gain > maxgain) {
        				maxgain = gain;
        				maxid = id;
        			}
        			
         		}
        		
        		if(maxgain <0 || maxid < 0) return ;
        		
        		selectPatternIndex.add(maxid);
        		
        		coverededges_curall.addAll(CoveredEdges_patterns.get(maxid));
        		//System.out.println(maxid + " is added, coverededges_curall.size():" + coverededges_curall.size());
        		
        		tempcount++;
        	}
        	
        	 int countofpatterns = 0;
             for(long i: selectPatternIndex) {
                Graph g  = allGraphs.get((int) i);
             	newreport(g, countofpatterns++);
             }
             int countofcoverededges = 0;
             Set<Integer> temp = new HashSet<Integer>();
             for(long i: selectPatternIndex) {
            	temp.addAll(CoveredEdges_patterns.get(i));
             }
             countofcoverededges = temp.size();
             
             System.out.println("Greedy,  Number of covered edges: " + countofcoverededges);
             int totalegdes = 0;
             for(int i=0;i< TRANS.size();i++) {
             	totalegdes += TRANS.get(i).getEdgeSize();
              }
             System.out.println("totalegdes : " +  totalegdes);
             System.out.println("Coverage rate : " + countofcoverededges*1.0 / totalegdes);
        }
    }
    private void  reportIndexSize(){
    	File outFile = new File("IndexTest");
        try (FileWriter writer = new FileWriter(outFile)) {
        	 for (Map.Entry<Integer, Set<Integer>> entry : Rpriv_i.entrySet()) {
        		 for(Integer e: entry.getValue())
        		    writer.write(e);
        		 //writer.write("\r\n");
        	 }
        	// writer.write("***********\r\n");
        	 writer.write(Priv_pattern.toString()+"\r\n");
        	 for (Map.Entry<Integer, Set<Integer>> entry : Rcov_edge.entrySet()) {
        		 for(Integer e: entry.getValue())
         		    writer.write(e);
         		// writer.write("\r\n");
        	 }
        	// writer.write("***********\r\n");
        	
           // private HashMap<Integer, Set<Integer>> CoveredEdges_OriginalGraphs;
        	 for (Map.Entry<Integer, Set<Integer>> entry : CoveredEdges_OriginalGraphs.entrySet()) {
        		 //writer.write("ssss" + "\r\n");
        		 for(Integer e: entry.getValue())
         		    writer.write(e);
         		// writer.write("\r\n");
        	 }
        	 
            //private Set<Integer>            allCoveredEdges;
            //private ArrayList<Set<Integer>> CoveredEdges_EachPattern;
            //// |priv(pattern)|: private edges for each pattern
           
            	
        } catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
        
      //  System.out.println("Index Maintenance Time: " + arg.maintainTime);
    }
    
    public int getLossScore(Set<Integer> dropededges, Long deleteid) {
    	int loss_count = 0;
    	
    	Set<Integer> set_temp  =  new HashSet<Integer>();
    	for (Long key : CoveredEdges_patterns.keySet()) {
    		 if(key != deleteid)   
    			 set_temp.addAll(CoveredEdges_patterns.get(key));
 		}
    	Iterator it = dropededges.iterator();
        while(it.hasNext()){
        	Integer edgeid  = (Integer)it.next();
        	if(!set_temp.contains(edgeid))  {
        		loss_count++;
        	}
        }
		return loss_count;
    }
    
    private void read(FileReader is) throws IOException {
    	int count = 0;
        BufferedReader read = new BufferedReader(is);
        while (true) {
            Graph g = new Graph(directed);
            read = g.read(read);
            if (g.isEmpty())
                break;
            TRANS.add(g);
            count++;
            if(count == arg.numberofgraphs) break; 
        }
        read.close();
    }
    
    private void readLight(FileReader is) throws IOException {
    	ArrayList<Graph> TRANS_Batch = new ArrayList<Graph>();
    	int iter = 0;
    	int count = 0;
        BufferedReader read = new BufferedReader(is);
        while (true) {
            Graph g = new Graph(directed);
            read = g.read(read);
            if (g.isEmpty())
                break;
            
            if(count < arg.ReadNumInEachBatch*(iter+1)) {
        		TRANS_Batch.add(g);
        		if(count == arg.numberofgraphs -1) {
            		BatchProcessor(TRANS_Batch);
            	}
        	}
        	else if(count == arg.ReadNumInEachBatch*(iter+1)) {
        		BatchProcessor(TRANS_Batch);
        		iter++;
        		TRANS_Batch.clear();
        		TRANS_Batch.add(g);
        	}
            count++;
            if(count == arg.numberofgraphs) break; 
        }
        read.close();
    }
    private void readLight2(FileReader is) throws IOException {
    	ArrayList<Graph> TRANS_Temp = new ArrayList<Graph>();
    	double MaxE, MaxV, AvgE, AvgV;
    	MaxE = MaxV = AvgE = AvgV = 0;
    	int count = 0;
        BufferedReader read = new BufferedReader(is);
        while (true) {
            Graph g = new Graph(directed);
            read = g.read(read);
            if (g.isEmpty())
                break;
            //TRANS.add(g);
            TRANS_Temp.add(g);
            AvgE += g.getEdgeSize();
            AvgV += g.size();
            if(g.getEdgeSize() > MaxE ) MaxE = g.getEdgeSize();
            if(g.size()  > MaxV )       MaxV = g.size();
            
            count++;
            if(count == arg.numberofgraphs) break; 
        }
        read.close();
        System.out.println(count + " graphs are readed in." ) ;
        System.out.println(AvgE + " edges are readed in." ) ;
        System.out.println(AvgV + " nodes are readed in." ) ;
        System.out.println(MaxE + "," + MaxV + "," + AvgE/count + "," +  AvgV/count) ;
        
        ArrayList<Graph> TRANS_Batch = new ArrayList<Graph>();
        int iter = 0;
        for(int i= 0; i< TRANS_Temp.size();i++) {
        	if(i < arg.ReadNumInEachBatch*(iter+1)) {
        		TRANS_Batch.add(TRANS_Temp.get(i));
        		if(i == TRANS_Temp.size() -1) {
            		BatchProcessor(TRANS_Batch);
            	}
        	}
        	else if(i == arg.ReadNumInEachBatch*(iter+1)) {
        		BatchProcessor(TRANS_Batch);
        		iter++;
        		TRANS_Batch.clear();
        		TRANS_Batch.add(TRANS_Temp.get(i));
        	}
        }
        
    }
    private void  BatchProcessor( ArrayList<Graph> TRANS_Batch) throws IOException {
    	System.out.println("TRANS_Batch:" + TRANS_Batch.size());
    	// Selected top VertorSize edges to represent a graph
    	int VertorSize = 10;
    	ArrayList<ArrayList<Double>> featurevectors = new ArrayList<ArrayList<Double>>();
    	for(int i=0;i<TRANS_Batch.size();i++) {
    		ArrayList<Double> temp = new ArrayList<Double>();
    		for(int j=0;j<VertorSize;j++) {
    			temp.add(0.0);
    		}
    		featurevectors.add(temp);
    	}
    	
    	
    	//step 1: extract all edges 
    	
    	//singleEdge[graphid][edge][count]
    	//singleEdgeLabel[edge][count]
    	NavigableMap<Integer, NavigableMap<SimpleEdge, Integer>> singleEdge= new TreeMap<>();
        NavigableMap<SimpleEdge, Integer> singleEdgeLabel= new TreeMap<>();
        for (int id = 0; id < TRANS_Batch.size(); ++id) {
            Graph g = TRANS_Batch.get(id);
            for (int from = 0; from < g.size(); ++from) {
           	 ArrayList<Edge> edges = new ArrayList<>();
                if (Misc.getForwardRoot(g, g.get(from), edges)) {
                    for (Edge it : edges) {
                        int key_1 = g.get(from).label;
                        int key_2 = it.eLabel;
                        int key_3 = g.get(it.to).label;
                        SimpleEdge e = new SimpleEdge();
                        e.vlabel1 = key_1;
                        e.vlabel2 = key_3;
                        if(key_1 == key_3 && it.from > it.to) continue; 
                        singleEdge.computeIfAbsent(id, k -> new TreeMap<>());
                        if(singleEdge.get(id).containsKey(e) == false) {
                        	singleEdgeLabel.put(e, Common.getValue(singleEdgeLabel.get(e))+1);
                        }
                        singleEdge.get(id).put(e, Common.getValue(singleEdge.get(id).get(e)) + 1);
                    }
                }
            }
        }
        List<Map.Entry<SimpleEdge, Integer>> entryArrayList = new ArrayList<>(singleEdgeLabel.entrySet());
        Collections.sort(entryArrayList, Comparator.comparing(Map.Entry::getValue));
       
        //step 2: select top VertorSize edges  
        ArrayList<SimpleEdge> singleEdgeLabel_selected = new ArrayList<>();
        
        if(entryArrayList.size()<VertorSize) VertorSize = entryArrayList.size();
        int count = 0;
        for (Map.Entry<SimpleEdge, Integer> entry : entryArrayList) {
        	count++;
        	if(entryArrayList.size() < (VertorSize + count)) {
        		singleEdgeLabel_selected.add(entry.getKey());
        		//System.out.println(entry.getKey().vlabel1 + "," + entry.getKey().vlabel2 + " freq: " + entry.getValue() );
        	}
        }
        
        //step 3: represent each graph as a feature vector
        for(Map.Entry<Integer, NavigableMap<SimpleEdge, Integer>> entry:singleEdge.entrySet())
        {
            int gid = entry.getKey();
            int edgesize = TRANS_Batch.get(gid).getEdgeSize();
            //System.out.println("***gid:" + gid + " , egdesize:"+edgesize);
        	NavigableMap<SimpleEdge, Integer> covered = entry.getValue();
        	//for (Map.Entry<SimpleEdge, Integer> temp : covered.entrySet()) {
        	//	 System.out.println(temp.getKey().vlabel1 + "," + temp.getKey().vlabel2 + " freq: " + temp.getValue());
           // }
        	int i = 0;
        	for(SimpleEdge e : singleEdgeLabel_selected) {
        		int freq = 0;
        		if(covered.get(e) != null) {
        			freq = covered.get(e);
        		}
        		featurevectors.get(gid).set(i, freq*1.0/edgesize);
        		i++;
        	}
        }
        //System.out.println("***Print featurevectors***");
        //for(int i=0;i<featurevectors.size();i++) {
        //	for(int j=0;j<featurevectors.get(i).size();j++) {
        //		System.out.print(featurevectors.get(i).get(j) + " ");
        //	}
        //	System.out.println();
        //}
        
       //step 4: clustering 
        int numPoints = TRANS_Batch.size();
        int dimensions = VertorSize;
        int numberofcluster = arg.numberofpatterns;
        double[][] points = new double[numPoints][dimensions];
        for(int i=0;i<numPoints;i++) {
        	for(int j=0;j<VertorSize;j++) {
        		points[i][j] = featurevectors.get(i).get(j);
        	}
        }
        
        final long startTime = System.currentTimeMillis();
        KMeans clustering = new KMeans.Builder(numberofcluster, points).iterations(1).pp(true).epsilon(.0001).useEpsilon(true).build();
        final long endTime = System.currentTimeMillis();
        // print timing information
        final long elapsed = endTime - startTime;
        System.out.println("Clustering took " + (double) elapsed/1000 + " seconds");
        // get output
        double[][] centroids = clustering.getCentroids();
        double WCSS          = clustering.getWCSS();
        int[]  assignment    = clustering.getAssignment();    
       // System.out.println("The centroids is: ");
       // for (int i = 0; i < centroids.length; i++) {
       // 	for(int j=0;j<centroids[0].length;j++) {
       // 		System.out.print(centroids[i][j] + "  ");
       // 	}
       // 	System.out.println();
       // }
      //  System.out.println("The within-cluster sum-of-squares (WCSS) = " + WCSS);
       // System.out.println("The assignment is: ");
       // for (int i = 0; i < assignment.length; i++) {
       // 	System.out.println("graph "+ i + " is in cluster "+assignment[i]);
       // }
        ArrayList<ArrayList<Integer>> clusters = new  ArrayList<ArrayList<Integer>>();
        for(int i= 0; i< numberofcluster;i++) {
        	ArrayList<Integer> temp = new ArrayList<Integer>();
        	for(int j=0;j< assignment.length;j++) {
        		if(i==assignment[j]) {
        			temp.add(j);
        		}
        	}
        	clusters.add(temp);
        }
        for(int i=0;i<clusters.size();i++) {
        	ArrayList<Integer> temp = clusters.get(i);
        	System.out.println("temp.size():" + temp.size());
        	int tempsize = (int)(temp.size() * (arg.SampleEdgeNum / arg.AvgE) * 1.0 / arg.numberofgraphs ); 
        	System.out.println("tempsize:" + tempsize);
        	Collections.shuffle(temp);
        	for(int k=0;k<tempsize;k++) {
        		this.TRANS.add(TRANS_Batch.get(temp.get(k)));
        	}
        }
       
    }
    
    private void  InitialPatternGenerator() throws IOException {
        	 ArrayList<Edge> edges = new ArrayList<>();
             NavigableMap<Integer, NavigableMap<Integer, NavigableMap<Integer, Projected>>> root = new TreeMap<>();
             for (int id = 0; id < TRANS.size(); ++id) {
                 Graph g = TRANS.get(id);
                 for (int from = 0; from < g.size(); ++from) {
                     if (Misc.getForwardRoot(g, g.get(from), edges)) {
                         for (Edge it : edges) {
                             int key_1 = g.get(from).label;
                             NavigableMap<Integer, NavigableMap<Integer, Projected>> root_1 = root.computeIfAbsent(key_1, k -> new TreeMap<>());
                             int key_2 = it.eLabel;
                             NavigableMap<Integer, Projected> root_2 = root_1.computeIfAbsent(key_2, k -> new TreeMap<>());
                             int key_3 = g.get(it.to).label;
                             Projected root_3 = root_2.get(key_3);
                             if (root_3 == null) {
                                 root_3 = new Projected();
                                 root_2.put(key_3, root_3);
                             }
                             root_3.push(id, it, null);
                         }
                     }
                 }
             }
             for (Entry<Integer, NavigableMap<Integer, NavigableMap<Integer, Projected>>> fromLabel : root.entrySet()) {
                 for (Entry<Integer, NavigableMap<Integer, Projected>> eLabel : fromLabel.getValue().entrySet()) {
                     for (Entry<Integer, Projected> toLabel : eLabel.getValue().entrySet()) {
                         DFS_CODE.push(0, 1, fromLabel.getKey(), eLabel.getKey(), toLabel.getKey());
                         project_Initial(toLabel.getValue());
                         DFS_CODE.pop();
                     }
                 }
             }
             System.out.println("After Initial Swapping, Number of covered edges: " + allCoveredEdges.size());
           	 int totalegdes = 0;
             for(int i=0;i< TRANS.size();i++) {
                totalegdes += TRANS.get(i).getEdgeSize();
             }
             System.out.println("totalegdes : " +  totalegdes);
             System.out.println("Coverage rate : " + allCoveredEdges.size()*1.0 / totalegdes);
    }
   
    private void  runIntern() throws IOException {
        // In case 1 node sub-graphs should also be mined for, do this as pre-processing step.
        if ( arg.minNodeNum <= 1) {
            /*
             * Do single node handling, as the normal gSpan DFS code based
             * processing cannot find sub-graphs of size |sub-g|==1. Hence, we
             * find frequent node labels explicitly.
             */
            for (int id = 0; id < TRANS.size(); ++id) {
                for (int nid = 0; nid < TRANS.get(id).size(); ++nid) {
                    int key = TRANS.get(id).get(nid).label;
                    //note: if singleVertex.get(id) is null, assign new TreeMap<>() to the key, ie, id
                    singleVertex.computeIfAbsent(id, k -> new TreeMap<>());
                    if (singleVertex.get(id).get(key) == null) {
                        // number of graphs it appears in
                        singleVertexLabel.put(key, Common.getValue(singleVertexLabel.get(key)) + 1);
                    }
                    singleVertex.get(id).put(key, Common.getValue(singleVertex.get(id).get(key)) + 1);
                }
            }
        }
        /*
         * All minimum support node labels are frequent 'sub-graphs'.
         * singleVertexLabel[nodeLabel] gives the number of graphs it appears in.
         */
        for (Entry<Integer, Integer> it : singleVertexLabel.entrySet()) {
            if (it.getValue() < arg.minSup)
                continue;

            int frequent_label = it.getKey();

            // Found a frequent node label, report it.
            Graph g = new Graph(directed);
            Vertex v = new Vertex();
            v.label = frequent_label;
            g.add(v);

            // [graph_id] = count for current substructure
            Vector<Integer> counts = new Vector<>();
            counts.setSize(TRANS.size());
            for (Entry<Integer, NavigableMap<Integer, Integer>> it2 : singleVertex.entrySet()) {
                counts.set(it2.getKey(), it2.getValue().get(frequent_label));
            }

            NavigableMap<Integer, Integer> gyCounts = new TreeMap<>();
            for (int n = 0; n < counts.size(); ++n)
                gyCounts.put(n, counts.get(n));

            reportSingle(g, gyCounts);
        }

        ArrayList<Edge> edges = new ArrayList<>();
        // note: [vertex1.label][eLabel][vertex2.label] = Projected
        NavigableMap<Integer, NavigableMap<Integer, NavigableMap<Integer, Projected>>> root = new TreeMap<>();

        for (int id = 0; id < TRANS.size(); ++id) {
            Graph g = TRANS.get(id);
            for (int from = 0; from < g.size(); ++from) {
                if (Misc.getForwardRoot(g, g.get(from), edges)) {
                    for (Edge it : edges) {
                        int key_1 = g.get(from).label;
                        NavigableMap<Integer, NavigableMap<Integer, Projected>> root_1 = root.computeIfAbsent(key_1, k -> new TreeMap<>());
                        int key_2 = it.eLabel;
                        NavigableMap<Integer, Projected> root_2 = root_1.computeIfAbsent(key_2, k -> new TreeMap<>());
                        int key_3 = g.get(it.to).label;
                        Projected root_3 = root_2.get(key_3);
                        if (root_3 == null) {
                            root_3 = new Projected();
                            root_2.put(key_3, root_3);
                        }
                        root_3.push(id, it, null);
                    }
                }
            }
        }

        for (Entry<Integer, NavigableMap<Integer, NavigableMap<Integer, Projected>>> fromLabel : root.entrySet()) {
            for (Entry<Integer, NavigableMap<Integer, Projected>> eLabel : fromLabel.getValue().entrySet()) {
                for (Entry<Integer, Projected> toLabel : eLabel.getValue().entrySet()) {
                    // Build the initial two-node graph. It will be grown recursively within project.
                	//note: 0,1, vertex1_label, eLabel, vertex2_label
                    DFS_CODE.push(0, 1, fromLabel.getKey(), eLabel.getKey(), toLabel.getKey());
                   //note: The position of  edge(vertex1_label, eLabel, vertex2_label) occurs in.  
                   //It contains Projected object with a set of PDFS where each of them contains an edge in original graph < graph id, the original edge in the graph,  prev = null>
                    project(toLabel.getValue());
                    //System.out.println("allCoveredEdges:" + this.allCoveredEdges.size()+"," + this.allGraphs.size());
                    DFS_CODE.pop();
                }
            }
        }
    }

    private void reportSingle(Graph g, NavigableMap<Integer, Integer> nCount) throws IOException {
        int sup = 0;
        
        // note: total occurrences, [graph] = nCount
        for (Entry<Integer, Integer> it : nCount.entrySet()) {
            sup += Common.getValue(it.getValue());
        }

        if ( arg.maxNodeNum  > arg.minNodeNum && g.size() > arg.maxNodeNum)
            return;
        if (arg.minNodeNum > 0 && g.size() < arg.minNodeNum)
            return;

        os.write("t # " + ID + " * " + sup + System.getProperty("line.separator"));
        g.write(os);
        ID++;
    }

    
    private void newreport(Graph g, int id) throws IOException {
        // Filter to small/too large graphs.
    	int sup = 0;
    	int ID = id;
      //  Graph g = new Graph(directed);
      //  code.toGraph(g);
        os.write("Final t # " + ID + " * " + sup + System.getProperty("line.separator"));
        g.write(os);
        ++ID;
    }
    private void newreportbeforeswap(Graph g, int id) throws IOException {
        // Filter to small/too large graphs.
    	int sup = 0;
    	int ID = id;
      //  Graph g = new Graph(directed);
      //  code.toGraph(g);
        os.write("BeforeSwap t # " + ID + " * " + sup + System.getProperty("line.separator"));
        g.write(os);
        ++ID;
    }
    
   
    public int getBenefitScore(Set<Integer> coverededges) {
    	int unique_count = 0;
    	Iterator it = coverededges.iterator();
        while(it.hasNext()){
        	if(allCoveredEdges.contains(it.next())) continue;
        	unique_count++;
        }
		return unique_count;
    }
    
    
    public void Delete(int deleteid) {
    	//Long Time1 = System.currentTimeMillis();
    	//if(numberofcovered != allCoveredEdges.size())
    	//{System.out.println("error!!!!!!!!!!!!"); return;}
    	//System.out.println("##########delete##############");
    	
    	//allGraphs.remove(deleteid);
    	allGraphs.set(deleteid, null);
    	
    	Rpriv_i.get(Priv_pattern.get(deleteid)).remove(deleteid);
         
    	
    	
    	Set<Integer> coverededges_pattern = CoveredEdges_EachPattern.get(deleteid);
    	
    	for(Integer e: coverededges_pattern) {
    		
    		Rcov_edge.get(e).remove(deleteid);
    		
    		if(Rcov_edge.get(e).size() == 0) {
    			numberofcovered--;
    			
    			allCoveredEdges.remove(e);
    			
    			//CoveredEdges_OriginalGraphs.get(e/1000).remove(e);
    			
    		}else if (Rcov_edge.get(e).size() == 1) {
        		  int tempid = Rcov_edge.get(e).iterator().next();
        		  //System.out.println("tempid: " + tempid);
        		  int temp = Priv_pattern.get(tempid);
        		  Priv_pattern.set(tempid,  temp+1);
        		  Rpriv_i.get(temp).remove(tempid);
        		  if(Rpriv_i.get(temp+1) ==null) Rpriv_i.put(temp+1,  new HashSet<Integer>()); 
        		  Rpriv_i.get(temp+1).add(tempid);
    		}
    		
    		
    		
    	}
    	
    	
    	 
		// update  CoveredEdges_OriginalGraphs
		 CoveredEdges_OriginalGraphs.clear();
		 for(Integer edgeid : allCoveredEdges) {
			 Integer gid = edgeid / 1000;
			 Set<Integer>  temp = CoveredEdges_OriginalGraphs.get(gid);
			 if(temp == null) temp = new  HashSet<Integer>();
			 temp.add(edgeid);
			 CoveredEdges_OriginalGraphs.put(gid, temp);
		 }
		 
    	//Priv_pattern.remove(deleteid);
    	Priv_pattern.set(deleteid, -1);
    	
    	//CoveredEdges_EachPattern.remove(deleteid);
    	CoveredEdges_EachPattern.set(deleteid, new HashSet<Integer>());
    	
    	 // System.out.println("**************");
         // System.out.println("Rcov_edge: "+ Rcov_edge.toString());
         // System.out.println("Priv_pattern: "+  Priv_pattern.toString());
         // System.out.println("Rpriv_i: "+   Rpriv_i.toString());
        //  System.out.println("numberofcovered: "+ numberofcovered);
        //  System.out.println("minimumpattern_id: "+ minimumpattern_id);
        //  System.out.println("minimumpattern_score: "+ minimumpattern_score);
   // 	Long Time2 = System.currentTimeMillis();
	//	arg.maintainTime += (Time2 - Time1)*1.0 /1000;
    }
    
    
   public void Insert(Projected projected, int insertid) {
	   //Long Time1 = System.currentTimeMillis();
	   //System.out.println("Insert");
	    Graph g = new Graph(directed);
        DFS_CODE.toGraph(g);
        if(allGraphs.size() < arg.numberofpatterns) allGraphs.add(g);
        else allGraphs.set(insertid, g);
        Set<Integer> coverededges_pattern = new HashSet<Integer>();
        if(Priv_pattern.size() < arg.numberofpatterns) Priv_pattern.add(0);
        else Priv_pattern.set(insertid,  0);
        for (PDFS aProjected : projected) {
      	    int id = aProjected.id;
      	    //System.out.println("id:" + id);
      	    Set<Integer> tempedges = new HashSet<Integer>(); 
            for (PDFS p = aProjected; p != null; p = p.prev) {
          	  Integer temp = 1000 * id + p.edge.id;
          	  coverededges_pattern.add(temp);
          	  tempedges.add(temp);
          	 allCoveredEdges.add(temp);
            } 
            if(CoveredEdges_OriginalGraphs.containsKey(id)) {
         	   tempedges.addAll(CoveredEdges_OriginalGraphs.get(id));
               CoveredEdges_OriginalGraphs.replace(id, tempedges);
            }else {
         	   CoveredEdges_OriginalGraphs.put(id, tempedges);
             }
         }
        if(CoveredEdges_EachPattern.size() < arg.numberofpatterns) CoveredEdges_EachPattern.add(coverededges_pattern);
        else CoveredEdges_EachPattern.set(insertid, coverededges_pattern);
        for (int temp : coverededges_pattern ) {
        	if(Rcov_edge.get(temp) ==null) Rcov_edge.put(temp, new HashSet<Integer>()); 
           	Rcov_edge.get(temp).add(insertid);
           	if(Rcov_edge.get(temp).size() == 1) {
           		//edgenum++;
           		Priv_pattern.set(insertid,  Priv_pattern.get(insertid)+1);
           		numberofcovered++;
            }else if(Rcov_edge.get(temp).size() == 2) {
           		  int tempid = -1;
           		  for(int e : Rcov_edge.get(temp)) if(e != insertid)  tempid = e;
           		  Priv_pattern.set(tempid,  Priv_pattern.get(tempid)-1);
           		  
           		  Rpriv_i.get(Priv_pattern.get(tempid)+1).remove(tempid);
           		  if(Rpriv_i.get(Priv_pattern.get(tempid)) ==null) Rpriv_i.put(Priv_pattern.get(tempid),  new HashSet<Integer>()); 
           		  Rpriv_i.get(Priv_pattern.get(tempid)).add(tempid);
           	}
		}
       // System.out.println("%%%%%Priv_pattern: "+  Priv_pattern.toString());
        //System.out.println("%%%%%Rpriv_i: "+   Rpriv_i.toString());
        //Priv_pattern.add(edgenum);
        if(Rpriv_i.get(Priv_pattern.get(insertid)) ==null) Rpriv_i.put(Priv_pattern.get(insertid),  new HashSet<Integer>());  
        Rpriv_i.get(Priv_pattern.get(insertid)).add(insertid);
        for(int i=0;i<= Priv_pattern.get(insertid);i++) {
        	if(Rpriv_i.get(i) != null && Rpriv_i.get(i).size()>0) {
        		//System.out.println(Priv_pattern.get(insertid)+","+i+ "," + Rpriv_i.get(i).size());
        		minimumpattern_id = Rpriv_i.get(i).iterator().next();
        		minimumpattern_score =  Priv_pattern.get(minimumpattern_id);
        		break;
        	}
        }
       //System.out.println("############");
       // System.out.println("Rcov_edge: "+ Rcov_edge.toString());
       // System.out.println("Priv_pattern: "+  Priv_pattern.toString());
       // System.out.println("Rpriv_i: "+   Rpriv_i.toString());
       // System.out.println("numberofcovered: "+ numberofcovered);
       // System.out.println("minimumpattern_id: "+ minimumpattern_id);
       // System.out.println("minimumpattern_score: "+ minimumpattern_score);
        //++ID; 
       // System.out.println("allCoveredEdges.size(): " + allCoveredEdges.size());
        
        //Long Time2 = System.currentTimeMillis();
		//arg.maintainTime += (Time2 - Time1)*1.0 /1000;
    }
   
   public void InsertWithSimpleIndex(Projected projected, int insertid) {
   	   Graph g = new Graph(directed);
       DFS_CODE.toGraph(g);
       if(allGraphs.size() < arg.numberofpatterns) allGraphs.add(g);
       //else allGraphs.set(insertid, g);
       else if(allGraphs.size() > insertid) allGraphs.set(insertid, g);
       else allGraphs.add(g);
       Set<Integer> coverededges_pattern = new HashSet<Integer>();
       for (PDFS aProjected : projected) {
       	   int id = aProjected.id;
       	   Set<Integer> tempedges = new HashSet<Integer>(); 
             for (PDFS p = aProjected; p != null; p = p.prev) {
           	  Integer temp = 1000 * id + p.edge.id;
           	  coverededges_pattern.add(temp);
           	  if(PatternsID_edges.containsKey(temp)) {
           		 Set<Integer> temppatternids = PatternsID_edges.get(temp);
           		 temppatternids.add(id);
           		 PatternsID_edges.replace(temp, temppatternids);
           	  }else {
           		 Set<Integer> temppatternids =new HashSet<Integer>();
           		 temppatternids.add(id);	 
           		 PatternsID_edges.put(temp, temppatternids);
           	  }
           	 allCoveredEdges.add(temp);
           	 tempedges.add(temp);
             } 
             if(CoveredEdges_OriginalGraphs.containsKey(id)) {
          	   tempedges.addAll(CoveredEdges_OriginalGraphs.get(id));
                 CoveredEdges_OriginalGraphs.replace(id, tempedges);
             }else {
          	   CoveredEdges_OriginalGraphs.put(id, tempedges);
             }
               
       }
       CoveredEdges_patterns.put((long) insertid, coverededges_pattern);
       
       //System.out.println("allCoveredEdges Size: " + allCoveredEdges.size());
   }
    
    public  Boolean reportwithCoveredEdges(int sup, Projected projected) throws IOException {
    	// Filter to small/too large graphs.
        if (arg.maxNodeNum > arg.minNodeNum && DFS_CODE.countNode() > arg.maxNodeNum)
            return false;
        if (arg.minNodeNum > 0 && DFS_CODE.countNode() < arg.minNodeNum )
            return false;
        //////////////////////////////////
        if( arg.strategy.equals("greedy") || allGraphs.size() <  arg.numberofpatterns) {
        	
        	 if(arg.isPESIndex) Insert(projected, allGraphs.size());
        	 else if(arg.isSimpleIndex)  InsertWithSimpleIndex(projected, allGraphs.size()) ;
        	 else ;
        
             if(allGraphs.size()  == arg.numberofpatterns) {
            	  int count = 0;
                  for(Graph tempg : allGraphs) {
                  	newreportbeforeswap(tempg, count++);
                  }
            	 System.out.println("Before Swapping, Number of covered edges: " + allCoveredEdges.size());
            	 int totalegdes = 0;
                 for(int i=0;i< TRANS.size();i++) {
                 	totalegdes += TRANS.get(i).getEdgeSize();
                  }
                 System.out.println("totalegdes : " +  totalegdes);
                 System.out.println("Coverage rate : " + allCoveredEdges.size()*1.0 / totalegdes);
                 
                 if(arg.isSimpleIndex && arg.strategy.equals("topk")) {
                	 int patternid_min = 0;
                	 int loss_score_min = Integer.MAX_VALUE;
                	 for (Long key : CoveredEdges_patterns.keySet()) {
                		 Set<Integer> dropededges =  CoveredEdges_patterns.get(key);
                		 int loss_score = getLossScore(dropededges, key);
                		 if(loss_score < loss_score_min) {
                			 loss_score_min  =  loss_score;
                			 patternid_min = key.intValue();
                		 } 
                	 }
                	 minimumpattern_score  =  loss_score_min;
                	 minimumpattern_id = patternid_min;
                	 
                	 //System.out.println("minimumpattern_score: " + minimumpattern_score);
                	 //System.out.println("minimumpattern_id: " + minimumpattern_id);
                 }
             }
        }
        else {
        	 /// 1. calculate benefit score
        	 int benefit_score = 0;
        	 //// 2. calculate minimum loss score
        	 int patternid_min = -1;
        	 int loss_score_min = Integer.MAX_VALUE;
        	 
        	 if(arg.isPESIndex) {
        		 Set<Integer> coverededges_pattern = new HashSet<Integer>();
                 for (PDFS aProjected : projected) {
                 	  int id = aProjected.id;
                       for (PDFS p = aProjected; p != null; p = p.prev) {
                     	  Integer temp = 1000 * id + p.edge.id;
                     	  coverededges_pattern.add(temp);
                       } 
                 }
                 for(Integer e : coverededges_pattern) {
                	 if( this.Rcov_edge.get(e) == null || this.Rcov_edge.get(e).size() == 0) {
                		 benefit_score++;
                	 }
                 }
                 patternid_min = minimumpattern_id;
            	 loss_score_min = minimumpattern_score;
            	 
            	 Boolean swapflag = false;
            	 if(arg.swapcondition.equals("swap1")) {
            		 if(benefit_score > 2* loss_score_min ) {
            			 swapflag = true;
                	 }
            	 }else  if(arg.swapcondition.equals("swap2")) {
            		 if(benefit_score > loss_score_min  + this.allCoveredEdges.size()*1.0/arg.numberofpatterns) {
            			 swapflag = true;
                	 }
            	 }else {
            		 if(benefit_score > (1+arg.swapAlpha)*loss_score_min  + (1-arg.swapAlpha)*(this.allCoveredEdges.size()*1.0/arg.numberofpatterns)) {
            			 swapflag = true;
                	 }
            	 }
            	 if(swapflag) {
            		 os.write(patternid_min + " is swapped out!");
            		 os.write("(Swapping Phase) benefit_score: " + benefit_score + ", loss_score_min: " + loss_score_min + "\\n");
            		 allGraphs.get(patternid_min).write(os);
            		 
            		 Delete(patternid_min);
            		 Insert(projected, patternid_min);
            		 
            		 
            	 }
            	 return swapflag;
        	 }else if(arg.isSimpleIndex) {
        		 Set<Integer> coverededges_pattern = new HashSet<Integer>();
                 for (PDFS aProjected : projected) {
                 	  int id = aProjected.id;
                       for (PDFS p = aProjected; p != null; p = p.prev) {
                     	  Integer temp = 1000 * id + p.edge.id;
                     	  coverededges_pattern.add(temp);
                       } 
                 }
            	 benefit_score = getBenefitScore(coverededges_pattern);
            	 
            	 loss_score_min  =  minimumpattern_score;
        		 patternid_min   =  minimumpattern_id ;
        		 
        		 Boolean swapflag = false;
            	 if(arg.swapcondition.equals("swap1")) {
            		 if(benefit_score > 2* loss_score_min ) {
            			 swapflag = true;
                	 }
            	 }else  if(arg.swapcondition.equals("swap2")) {
            		 if(benefit_score > loss_score_min  + this.allCoveredEdges.size()*1.0/arg.numberofpatterns) {
            			 swapflag = true;
                	 }
            	 }else {
            		 if(benefit_score > (1+arg.swapAlpha)*loss_score_min  + (1-arg.swapAlpha)*(this.allCoveredEdges.size()*1.0/arg.numberofpatterns)) {
            			 swapflag = true;
                	 }
            	 }
            	 if(swapflag) {
            		 os.write(patternid_min + " is swapped out!");
            		 os.write("(Swapping Phase) benefit_score: " + benefit_score + ", loss_score_min: " + loss_score_min );
            		 allGraphs.get(patternid_min).write(os);
            		 
            		 // update allAllGraphs 
            		 Graph g = new Graph(directed);
                     DFS_CODE.toGraph(g);
            		 allGraphs.set(patternid_min, g);
            		 
            
            		 // update CoveredEdges_patterns 
            		 CoveredEdges_patterns.replace((long) patternid_min, coverededges_pattern);
            		 
            		 // update allCoveredEdges 
            		 allCoveredEdges.clear();
            		 for (Long key : CoveredEdges_patterns.keySet()) {
            			 Set<Integer> temp = CoveredEdges_patterns.get(key);
            			 allCoveredEdges.addAll(temp);
            		 }
            		 
            		// update  CoveredEdges_OriginalGraphs
            		 CoveredEdges_OriginalGraphs.clear();
            		 for(Integer edgeid : allCoveredEdges) {
            			 Integer gid = edgeid / 1000;
            			 Set<Integer>  temp = CoveredEdges_OriginalGraphs.get(gid);
            			 if(temp == null) temp = new  HashSet<Integer>();
            			 temp.add(edgeid);
            			 CoveredEdges_OriginalGraphs.put(gid, temp);
            		 }
            		 
            		 
            		 // update PatternsID_edges 
            		 PatternsID_edges.clear();
            		 for(Integer edgeid : allCoveredEdges) {
            			 Set<Integer> tempedges =  new HashSet<Integer>();
            			 for (Long key : CoveredEdges_patterns.keySet()) {
                			 Set<Integer> temp = CoveredEdges_patterns.get(key);
                			 if(temp.contains(edgeid)) {
                				 tempedges.add(key.intValue());
                			 }
                		 }
            			 PatternsID_edges.put(edgeid, tempedges);
            		 }
            		 
            		 //System.out.println("allCoveredEdges Size: " + allCoveredEdges.size());
            		 if(true) {
                         patternid_min = 0;
                    	 loss_score_min = Integer.MAX_VALUE;
                    	 for (Long key : CoveredEdges_patterns.keySet()) {
                    		 Set<Integer> dropededges =  CoveredEdges_patterns.get(key);
                    		 int loss_score =  getLossScore(dropededges, key);
                    		 if(loss_score < loss_score_min) {
                    			 loss_score_min  =  loss_score;
                    			 patternid_min = key.intValue();
                    		 } 
                    	 }
                    	 minimumpattern_score  =  loss_score_min;
                    	 minimumpattern_id = patternid_min;
                    	 
                    	 //System.out.println("minimumpattern_score: " + minimumpattern_score);
                    	 //System.out.println("minimumpattern_id: " + minimumpattern_id);
                     }
            	 }
            	 return swapflag;
        	 }else {
        		 
        	 }
        }
        return false;
    }
    
    public  void reportwithCoveredEdges_Initial(int sup, Projected projected) throws IOException {
    	// Filter to small/too large graphs.
        if (arg.maxNodeNum > arg.minNodeNum && DFS_CODE.countNode() > arg.maxNodeNum)
            return;
        if (arg.minNodeNum > 0 && DFS_CODE.countNode() < arg.minNodeNum)
            return;
        
        //////////////////////////////////
        if(allGraphs.size() <  arg.numberofpatterns) {
        	
      	  Graph g = new Graph(directed);
          DFS_CODE.toGraph(g);
          os.write("Initial****t # " + allGraphs.size() + " * " + sup + System.getProperty("line.separator"));
          g.write(os);
            
       	 if(arg.isPESIndex) Insert(projected, allGraphs.size());
    	 else if(arg.isSimpleIndex)  InsertWithSimpleIndex(projected, allGraphs.size()) ;
    	 else ;
       	 //System.out.println("allCoveredEdges.size(): " + allCoveredEdges.size());
         if(allGraphs.size()  == arg.numberofpatterns) {
        	  int count = 0;
              for(Graph tempg : allGraphs) {
              	newreportbeforeswap(tempg, count++);
              }
        	 System.out.println("Before Swapping, Number of covered edges: " + allCoveredEdges.size());
        	 int totalegdes = 0;
             for(int i=0;i< TRANS.size();i++) {
             	totalegdes += TRANS.get(i).getEdgeSize();
              }
             System.out.println("totalegdes : " +  totalegdes);
             System.out.println("Coverage rate : " + allCoveredEdges.size()*1.0 / totalegdes);
             
             if(arg.isSimpleIndex && arg.strategy.equals("topk")) {
            	 int patternid_min = 0;
            	 int loss_score_min = Integer.MAX_VALUE;
            	 for (Long key : CoveredEdges_patterns.keySet()) {
            		 Set<Integer> dropededges =  CoveredEdges_patterns.get(key);
            		 int loss_score = getLossScore(dropededges, key);
            		 if(loss_score < loss_score_min) {
            			 loss_score_min  =  loss_score;
            			 patternid_min = key.intValue();
            		 } 
            	 }
            	 minimumpattern_score  =  loss_score_min;
            	 minimumpattern_id = patternid_min;
             }
         }
        }
        else {
       	 /// 1. calculate benefit score
       	 int benefit_score = 0;
       	 //// 2. calculate minimum loss score
       	 int patternid_min = -1;
       	 int loss_score_min = Integer.MAX_VALUE;
       	 
       	 if(arg.isPESIndex) {
       		 Set<Integer> coverededges_pattern = new HashSet<Integer>();
                for (PDFS aProjected : projected) {
                	  int id = aProjected.id;
                      for (PDFS p = aProjected; p != null; p = p.prev) {
                    	  Integer temp = 1000 * id + p.edge.id;
                    	  coverededges_pattern.add(temp);
                      } 
                }
                for(Integer e : coverededges_pattern) {
               	 if( this.Rcov_edge.get(e) == null || this.Rcov_edge.get(e).size() == 0) {
               		 benefit_score++;
               	 }
                }
                patternid_min = minimumpattern_id;
           	 loss_score_min = minimumpattern_score;
           	 
           	 Boolean swapflag = false;
           	 if(arg.swapcondition.equals("swap1")) {
           		 if(benefit_score > 2* loss_score_min ) {
           			 swapflag = true;
               	 }
           	 }else  if(arg.swapcondition.equals("swap2")) {
           		 if(benefit_score > loss_score_min  + this.allCoveredEdges.size()*1.0/arg.numberofpatterns) {
           			 swapflag = true;
               	 }
           	 }else {
           		 if(benefit_score > (1+arg.swapAlpha)*loss_score_min  + (1-arg.swapAlpha)*(this.allCoveredEdges.size()*1.0/arg.numberofpatterns)) {
           			 swapflag = true;
               	 }
           	 }
           	 if(swapflag) {
           		 os.write(patternid_min + " is swapped out!");
           		 os.write("Initial Swapping, benefit_score: " + benefit_score + ", loss_score_min: " + loss_score_min );
           		 allGraphs.get(patternid_min).write(os);
           		 
           		 Delete(patternid_min);
           		 Insert(projected, patternid_min);
           	 }
       	 }else if(arg.isSimpleIndex) {
       		 Set<Integer> coverededges_pattern = new HashSet<Integer>();
                for (PDFS aProjected : projected) {
                	  int id = aProjected.id;
                      for (PDFS p = aProjected; p != null; p = p.prev) {
                    	  Integer temp = 1000 * id + p.edge.id;
                    	  coverededges_pattern.add(temp);
                      } 
                }
           	 benefit_score = getBenefitScore(coverededges_pattern);
           	 
           	 loss_score_min  = minimumpattern_score ;
       		 patternid_min   = minimumpattern_id ;
       		 
       		 Boolean swapflag = false;
           	 if(arg.swapcondition.equals("swap1")) {
           		 if(benefit_score > 2* loss_score_min ) {
           			 swapflag = true;
               	 }
           	 }else  if(arg.swapcondition.equals("swap2")) {
           		 if(benefit_score > loss_score_min  + this.allCoveredEdges.size()*1.0/arg.numberofpatterns) {
           			 swapflag = true;
               	 }
           	 }else {
           		 if(benefit_score > (1+arg.swapAlpha)*loss_score_min  + (1-arg.swapAlpha)*(this.allCoveredEdges.size()*1.0/arg.numberofpatterns)) {
           			 swapflag = true;
               	 }
           	 }
           	 if(swapflag) {
           		 os.write(patternid_min + " is swapped out!");
           		 os.write("Initial Swapping, benefit_score: " + benefit_score + ", loss_score_min: " + loss_score_min );
           		 allGraphs.get(patternid_min).write(os);
           		 
           		 // update allAllGraphs 
           		 Graph g = new Graph(directed);
                    DFS_CODE.toGraph(g);
           		 allGraphs.set(patternid_min, g);
           		 
           
           		 // update CoveredEdges_patterns 
           		 CoveredEdges_patterns.replace((long) patternid_min, coverededges_pattern);
           		 
           		 // update allCoveredEdges 
           		 allCoveredEdges.clear();
           		 for (Long key : CoveredEdges_patterns.keySet()) {
           			 Set<Integer> temp = CoveredEdges_patterns.get(key);
           			 allCoveredEdges.addAll(temp);
           		 }
           		 
           		// update  CoveredEdges_OriginalGraphs
           		 CoveredEdges_OriginalGraphs.clear();
           		 for(Integer edgeid : allCoveredEdges) {
           			 Integer gid = edgeid / 1000;
           			 Set<Integer>  temp = CoveredEdges_OriginalGraphs.get(gid);
           			 if(temp == null) temp = new  HashSet<Integer>();
           			 temp.add(edgeid);
           			 CoveredEdges_OriginalGraphs.put(gid, temp);
           		 }
           		 
           		 
           		 // update PatternsID_edges 
           		 PatternsID_edges.clear();
           		 for(Integer edgeid : allCoveredEdges) {
           			 Set<Integer> tempedges =  new HashSet<Integer>();
           			 for (Long key : CoveredEdges_patterns.keySet()) {
               			 Set<Integer> temp = CoveredEdges_patterns.get(key);
               			 if(temp.contains(edgeid)) {
               				 tempedges.add(key.intValue());
               			 }
               		 }
           			 PatternsID_edges.put(edgeid, tempedges);
           		 }
           		 
           		 
           		 if(true) {
                     patternid_min = 0;
                   	 loss_score_min = Integer.MAX_VALUE;
                   	 for (Long key : CoveredEdges_patterns.keySet()) {
                   		 Set<Integer> dropededges =  CoveredEdges_patterns.get(key);
                   		 int loss_score =  getLossScore(dropededges, key);
                   		 if(loss_score < loss_score_min) {
                   			 loss_score_min  =  loss_score;
                   			 patternid_min = key.intValue();
                   		 } 
                   	 }
                   	 minimumpattern_score  =  loss_score_min;
                   	 minimumpattern_id = patternid_min;
                    }
           	 }
       	 }else {
       		 
       	 }	
       }
    }
    
    
    int DynamicSupportSetting() {
   	 //// calculate minimum loss score
   	 int loss_score_min =   this.minimumpattern_score;
     List<Integer> list = new ArrayList<Integer>();     
     for (int id = 0; id < TRANS.size(); ++id) {
    	 int count = 0;
         for (int nid = 0; nid < TRANS.get(id).size(); ++nid) {
        	for(Edge e : TRANS.get(id).get(nid).edge) {
        	    Integer edgeid = id*1000 + e.id;
        	    if(this.Rcov_edge.get(edgeid)==null || this.Rcov_edge.get(edgeid).size()==0) continue;
        	    count++;
        	}
         }
         list.add(count);
     }
     Collections.sort(list);
     int sum = 0;
     int index = list.size() -1;
     while(sum < 2 * loss_score_min) {
    	 sum += list.get(index);
    	 index--;
     }
     
     int ans = 1;
     
     if((list.size() - 1 - index) > ans) 
    	 ans = list.size() - 1 - index;
     if(ans <= (int)arg.minSup)  { ans = (int)arg.minSup; }
     //else {  System.out.println("DSS got better ans( > minSup): "+ ans) ; } 
     //System.out.println("DSS: "+ ans + ", loss_score_min: " + loss_score_min);
     return ans;
    }
    
    int DynamicSupportSetting2() {
      	 //// calculate minimum loss score
      	 int loss_score_min =   this.minimumpattern_score;
        List<Integer> list = new ArrayList<Integer>();     
        for (int id = 0; id < TRANS.size(); ++id) {
       	 int count = 0;
            for (int nid = 0; nid < TRANS.get(id).size(); ++nid) {
           	for(Edge e : TRANS.get(id).get(nid).edge) {
           	    Integer edgeid = id*1000 + e.id;
           	    if(this.allCoveredEdges.contains(edgeid)) continue;
           	    count++;
           	}
            }
            list.add(count);
        }
        Collections.sort(list);
        int sum = 0;
        int index = list.size() -1;
        while(sum < 2 * loss_score_min) {
       	 sum += list.get(index);
       	 index--;
        }
        
        int ans = 1;
        
        if((list.size() - 1 - index) > ans) 
       	 ans = list.size() - 1 - index;
        
        if(ans <= (int)arg.minSup)  {  ans = (int)arg.minSup;  }
       // else {  System.out.println("DSS got better ans( > minSup): "+ ans) ; } 
        
        
       // System.out.println("DSS: "+ ans + ", arg.minSup: " + arg.minSup);
        
        return ans;
      	 
      	 
       }

    Boolean BranchAndBound(Projected projected_g, Projected  projected_g2, Boolean hasupdated) {
    	//System.out.println("test");
    	//if(true) return false;
    	if(hasupdated) {
    		int maximum_benefit = 0;
    		Set<Integer> temp  = new HashSet<Integer>();
        	for (PDFS aProjected : projected_g2) {
             	   int id = aProjected.id;
             	   if(temp.contains(id)==false) {
             		  // E_i
             		  int size = this.TRANS.get(id).getEdgeSize();
             		  // Cov_i
             		  int count = 0;
             		  for(int i=0;i<size;i++) {
             			  Integer edgeid = 1000*id + i;
             			  //System.out.println("ss: " + id);
             			  if(CoveredEdges_OriginalGraphs.get(id) != null && CoveredEdges_OriginalGraphs.get(id).contains(edgeid)) count++;
             		  }
             		  maximum_benefit = maximum_benefit + size - count;
             		
             		 // if(maximum_benefit > 2*this.minimumpattern_score) {
                	//	   return false;
                	//   }
             		 if(arg.swapcondition.equals("swap1")) {
                		 if(maximum_benefit > 2* minimumpattern_score ) {
                			 return false;
                    	 }
                	 }else  if(arg.swapcondition.equals("swap2")) {
                		 if(maximum_benefit > minimumpattern_score  + this.allCoveredEdges.size()*1.0/arg.numberofpatterns) {
                			 return false;
                    	 }
                	 }else {
                		 if(maximum_benefit > (1+arg.swapAlpha)*minimumpattern_score  + (1-arg.swapAlpha)*(this.allCoveredEdges.size()*1.0/arg.numberofpatterns)) {
                			 return false;
                    	 }
                	 }
             		  
             		  temp.add(id);
             	   }
        	}
        	//System.out.println("test11111");
            return true;
    	}
    	else {
    		//System.out.println("test");
    		if(true) {
    			int maximum_benefit = 0;
        		Set<Integer> temp  = new HashSet<Integer>();
            	for (PDFS aProjected : projected_g2) {
                 	   int id = aProjected.id;
                 	   if(temp.contains(id)==false) {
                 		  // E_i
                 		  int size = this.TRANS.get(id).getEdgeSize();
                 		  // Cov_i
                 		  int count = 0;
                 		  if(CoveredEdges_OriginalGraphs.get(id) !=null) {
                 			 for(int i=0;i<size;i++) {
                    			  Integer edgeid = 1000*id + i;
                    			  if(CoveredEdges_OriginalGraphs.get(id).contains(edgeid)) count++;
                    		  }
                 		  }
                 		  maximum_benefit = maximum_benefit + size - count;
                 		  temp.add(id);
                 	   }
            	}
            	//if(maximum_benefit <= 2*this.minimumpattern_score)  {
            	//	return true;
            	//}
            	
            	 if(arg.swapcondition.equals("swap1")) {
            		 if(maximum_benefit <= 2* minimumpattern_score ) {
            			 return true;
                	 }
            	 }else  if(arg.swapcondition.equals("swap2")) {
            		 if(maximum_benefit <= minimumpattern_score  + this.allCoveredEdges.size()*1.0/arg.numberofpatterns) {
            			 return true;
                	 }
            	 }else {
            		 if(maximum_benefit <= (1+arg.swapAlpha)*minimumpattern_score  + (1-arg.swapAlpha)*(this.allCoveredEdges.size()*1.0/arg.numberofpatterns)) {
            			 return true;
                	 }
            	 }
    		}
    			
    			
    		int maximum_benefit = 0;
    		int totaledges = 0;
    		Set<Integer> graphIDs = new HashSet<Integer>();
        	Set<Integer> Cov_g  = new HashSet<Integer>();
        	Set<Integer> Cov_g2  = new HashSet<Integer>();
        	Set<Integer> Cov_i  = new HashSet<Integer>();
        	
        	//calculate Cov_g2
        	for (PDFS aProjected : projected_g2) {
          	    int id = aProjected.id;
                for (PDFS p = aProjected; p != null; p = p.prev) {
              	  Integer temp = 1000 * id + p.edge.id;
              	  Cov_g2.add(temp);
                } 
                graphIDs.add(id);
        	}
        	
        	//calculate Cov_g
        	for (PDFS aProjected : projected_g) {
          	    int id = aProjected.id;
          	    if(graphIDs.contains(id) == false) continue;
                for (PDFS p = aProjected; p != null; p = p.prev) {
              	  Integer temp = 1000 * id + p.edge.id;
              	  Cov_g.add(temp);
                } 
        	}
        	
        	//calculate Cov_i
        	for(int id: graphIDs) {
        		totaledges += this.TRANS.get(id).getEdgeSize();
        		//for(int e: allCoveredEdges) if(e >= 1000*id && e < 1000*(id+1)) Cov_i.add(e);
        		if(CoveredEdges_OriginalGraphs.get(id) == null) continue;
        		Cov_i.addAll(CoveredEdges_OriginalGraphs.get(id));
        	}
        	
        	 // Cov_diff = Cov(g) \ (Cov(g2) U Cov(g))
     	   Set<Integer> Cov_diff = new HashSet<Integer>();
     	   for(Integer e: Cov_g) {
     		   if(Cov_g2.contains(e) == false && Cov_i.contains(e) == false) {
     			   Cov_diff.add(e);
     		   }
     	   }
     	   
     	    // Cov_i U Cov_diff
     	   Set<Integer> Cov_union = Cov_i;
     	   Cov_union.addAll(Cov_diff);
     	   
     	   maximum_benefit =  totaledges - Cov_union.size();
     	   
 		 //  if(maximum_benefit > 2*this.minimumpattern_score) {
    	//	   return false;
    	 //  }
 		   
 		  if(arg.swapcondition.equals("swap1")) {
     		 if(maximum_benefit > 2* minimumpattern_score ) {
     			 return false;
         	 }
     	 }else  if(arg.swapcondition.equals("swap2")) {
     		 if(maximum_benefit > minimumpattern_score  + this.allCoveredEdges.size()*1.0/arg.numberofpatterns) {
     			 return false;
         	 }
     	 }else {
     		 if(maximum_benefit > (1+arg.swapAlpha)*minimumpattern_score  + (1-arg.swapAlpha)*(this.allCoveredEdges.size()*1.0/arg.numberofpatterns)) {
     			 return false;
         	 }
     	 }
 		   
 		   return true;
    	}
    }
    
    private void project(Projected projected) throws IOException {
   	  //Recursive sub-graph mining function (similar to sub-procedure 1 Sub-graph_Mining in [Yan2002]).
   	  //Check if the pattern is frequent enough.
    	int sup = support(projected);
       if (sup < arg.minSup)
           return;
       if (!isMin()) {
           return;
       }
       // Output the frequent substructure
       //  report(sup);
       Boolean hasupdated = reportwithCoveredEdges(sup, projected);
       //if(hasupdated) System.out.println("out");
       if(arg.maxNodeNum <= 2) return;
       
       
       //if(hasupdated) System.out.println("true"); else System.out.println("false");
       
       if(arg.hasDSS && hasupdated) {
       	// global pruning 
       	  arg.minSup = DynamicSupportSetting();
       }
       
       if (arg.maxNodeNum > arg.minNodeNum && DFS_CODE.countNode() > arg.maxNodeNum)
           return;

       /*
        * We just outputted a frequent sub-graph. As it is frequent enough, so
        * might be its (n+1)-extension-graphs, hence we enumerate them all.
        */
       ArrayList<Integer> rmPath = DFS_CODE.buildRMPath();
       int minLabel = DFS_CODE.get(0).fromLabel;
       int maxToc = DFS_CODE.get(rmPath.get(0)).to;

       
       //note:  1. new_bck_root[to_vertex][eLabel] =  Projected, since from_vertex is fixed as maxToc and labels of two end nodes are known. 
       //                 (Reason: for rightmost extension, backward edges can only be inserted when from_vertex is rightmost vertex)
       //            2. new_fwd_root[from_vertex][eLabel][to_vertex_label], since to_vertex is fixed as maxToc+1 and label of from_vertex is known. 
       //                 (Reason: for rightmost extension, forward edges can  be inserted with any nodes in the rightmost path)
       
       NavigableMap<Integer, NavigableMap<Integer, NavigableMap<Integer, Projected>>> new_fwd_root = new TreeMap<>();
       NavigableMap<Integer, NavigableMap<Integer, Projected>> new_bck_root = new TreeMap<>();
       ArrayList<Edge> edges = new ArrayList<>();

       // Enumerate all possible one edge extensions of the current substructure.
       for (PDFS aProjected : projected) {

           int id = aProjected.id;
           History history = new History(TRANS.get(id), aProjected);

           // XXX: do we have to change something here for directed edges?

           
           // backward
           for (int i = rmPath.size() - 1; i >= 1; --i) {
           	
           	//note: e1 = history.get(rmPath.get(i)),  e2 = history.get(rmPath.get(0)), check if there is an edge between e2.to and e1.from, and this edge is not already in history
               //           if yes, choose this edge as a backward edge
           	Edge e = Misc.getBackward(TRANS.get(id), history.get(rmPath.get(i)), history.get(rmPath.get(0)),history);
               if (e != null) {
                   int key_1 = DFS_CODE.get(rmPath.get(i)).from;
                   NavigableMap<Integer, Projected> root_1 = new_bck_root.computeIfAbsent(key_1, k -> new TreeMap<>());
                   int key_2 = e.eLabel;
                   Projected root_2 = root_1.get(key_2);
                   if (root_2 == null) {
                       root_2 = new Projected();
                       root_1.put(key_2, root_2);
                   }
                 //note:  the new Projected root_2 has a pointer to the previous Projected, aProjected
                   root_2.push(id, e, aProjected);
               }
           }

           // pure forward
           // FIXME: here we pass a too large e.to (== history[rmPath[0]].to
           // into getForwardPure, such that the assertion fails.
           //
           // The problem is:
           // history[rmPath[0]].to > TRANS[id].size()
           if (Misc.getForwardPure(TRANS.get(id), history.get(rmPath.get(0)), minLabel, history, edges))
               for (Edge it : edges) {
                   NavigableMap<Integer, NavigableMap<Integer, Projected>> root_1 = new_fwd_root.computeIfAbsent(maxToc, k -> new TreeMap<>());
                   int key_2 = it.eLabel;
                   NavigableMap<Integer, Projected> root_2 = root_1.computeIfAbsent(key_2, k -> new TreeMap<>());
                   int key_3 = TRANS.get(id).get(it.to).label;
                   Projected root_3 = root_2.get(key_3);
                   if (root_3 == null) {
                       root_3 = new Projected();
                       root_2.put(key_3, root_3);
                   }
                   root_3.push(id, it, aProjected);
               }
           // backtracked forward
           for (Integer aRmPath : rmPath)
               if (Misc.getForwardRmPath(TRANS.get(id), history.get(aRmPath), minLabel, history, edges))
                   for (Edge it : edges) {
                       int key_1 = DFS_CODE.get(aRmPath).from;
                       NavigableMap<Integer, NavigableMap<Integer, Projected>> root_1 = new_fwd_root.computeIfAbsent(key_1, k -> new TreeMap<>());
                       int key_2 = it.eLabel;
                       NavigableMap<Integer, Projected> root_2 = root_1.computeIfAbsent(key_2, k -> new TreeMap<>());
                       int key_3 = TRANS.get(id).get(it.to).label;
                       Projected root_3 = root_2.get(key_3);
                       if (root_3 == null) {
                           root_3 = new Projected();
                           root_2.put(key_3, root_3);
                       }
                       root_3.push(id, it, aProjected);
                   }
       }

       // Test all extended substructures.
       // backward
       for (Entry<Integer, NavigableMap<Integer, Projected>> to : new_bck_root.entrySet()) {
           for (Entry<Integer, Projected> eLabel : to.getValue().entrySet()) {
           	 if(arg.hasPRM  && this.allGraphs.size() == arg.numberofpatterns)
                {
                	if(BranchAndBound(projected, eLabel.getValue(),hasupdated)) {
                		continue;
                	}
                }
           	   DFS_CODE.push(maxToc, to.getKey(), -1, eLabel.getKey(), -1);
               project(eLabel.getValue());
               DFS_CODE.pop();
           }
       }

       //note:  There are many forward edges, so they should be visited descendingly resp.  from_vertex.
       //            e.g., let rightmost path be (1,2,4), then forward edges should be, (4,5), (2,5) and (1,5)
       // forward
       for (Entry<Integer, NavigableMap<Integer, NavigableMap<Integer, Projected>>> from : new_fwd_root.descendingMap().entrySet()) {
           for (Entry<Integer, NavigableMap<Integer, Projected>> eLabel : from.getValue().entrySet()) {
               for (Entry<Integer, Projected> toLabel : eLabel.getValue().entrySet()) {
               	if(arg.hasPRM  && this.allGraphs.size() == arg.numberofpatterns)
                   {
                   	if(BranchAndBound(projected, toLabel.getValue(),hasupdated)) {
                   		continue;
                   	}
                   }
               	   DFS_CODE.push(from.getKey(), maxToc + 1, -1, eLabel.getKey(), toLabel.getKey());
                   project(toLabel.getValue());
                   DFS_CODE.pop();
               }
           }
       }
   }
   
    private Integer getBenefitScore_Initial(Projected projected) {
    	//reportwithCoveredEdges_Initial(sup, projected);
        Set<Integer> coverededges_pattern = new HashSet<Integer>();
        for (PDFS aProjected : projected) {
        	  int id = aProjected.id;
              for (PDFS p = aProjected; p != null; p = p.prev) {
            	  Integer temp = 1000 * id + p.edge.id;
            	  coverededges_pattern.add(temp);
              } 
        }
   	    int benefitscore = getBenefitScore(coverededges_pattern);
   	    return benefitscore;
    }
    
    private void project_Initial(Projected projected) throws IOException {
        //int sup = support(projected);
       
    	//if(allDFSCodes.size() >=  numberofpatterns) {
    	//	return ;
    	//}
    	
        //if (!isMin()) {
        	//System.out.println("not minimum, number of nodes" +DFS_CODE.countNode()) ;
           // return;
        //}
        
        //reportwithCoveredEdges_Initial(sup, projected);
   	    int CurrentBenefitScore = getBenefitScore_Initial(projected);
        
   	    if(arg.maxNodeNum <= 2) return;
   	    
        if (arg.maxNodeNum > arg.minNodeNum && DFS_CODE.countNode() > arg.maxNodeNum)  return;
        ArrayList<Integer> rmPath = DFS_CODE.buildRMPath();
        int minLabel = DFS_CODE.get(0).fromLabel;
        int maxToc = DFS_CODE.get(rmPath.get(0)).to;
        NavigableMap<Integer, NavigableMap<Integer, NavigableMap<Integer, Projected>>> new_fwd_root = new TreeMap<>();
        NavigableMap<Integer, NavigableMap<Integer, Projected>> new_bck_root = new TreeMap<>();
        ArrayList<Edge> edges = new ArrayList<>();
        // Enumerate all possible one edge extensions of the current substructure.
        for (PDFS aProjected : projected) {
            int id = aProjected.id;
            History history = new History(TRANS.get(id), aProjected);
            // backward
            for (int i = rmPath.size() - 1; i >= 1; --i) {
            	//note: e1 = history.get(rmPath.get(i)),  e2 = history.get(rmPath.get(0)), check if there is an edge between e2.to and e1.from, and this edge is not already in history
                //           if yes, choose this edge as a backward edge
            	Edge e = Misc.getBackward(TRANS.get(id), history.get(rmPath.get(i)), history.get(rmPath.get(0)),history);
                if (e != null) {
                    int key_1 = DFS_CODE.get(rmPath.get(i)).from;
                    NavigableMap<Integer, Projected> root_1 = new_bck_root.computeIfAbsent(key_1, k -> new TreeMap<>());
                    int key_2 = e.eLabel;
                    Projected root_2 = root_1.get(key_2);
                    if (root_2 == null) {
                        root_2 = new Projected();
                        root_1.put(key_2, root_2);
                    }
                  //note:  the new Projected root_2 has a pointer to the previous Projected, aProjected
                    root_2.push(id, e, aProjected);
                }
            }
            // pure forward
            // FIXME: here we pass a too large e.to (== history[rmPath[0]].to
            // into getForwardPure, such that the assertion fails.
            //
            // The problem is:
            // history[rmPath[0]].to > TRANS[id].size()
            if (Misc.getForwardPure(TRANS.get(id), history.get(rmPath.get(0)), minLabel, history, edges))
                for (Edge it : edges) {
                    NavigableMap<Integer, NavigableMap<Integer, Projected>> root_1 = new_fwd_root.computeIfAbsent(maxToc, k -> new TreeMap<>());
                    int key_2 = it.eLabel;
                    NavigableMap<Integer, Projected> root_2 = root_1.computeIfAbsent(key_2, k -> new TreeMap<>());
                    int key_3 = TRANS.get(id).get(it.to).label;
                    Projected root_3 = root_2.get(key_3);
                    if (root_3 == null) {
                        root_3 = new Projected();
                        root_2.put(key_3, root_3);
                    }
                    root_3.push(id, it, aProjected);
                }
            // backtracked forward
            for (Integer aRmPath : rmPath)
                if (Misc.getForwardRmPath(TRANS.get(id), history.get(aRmPath), minLabel, history, edges))
                    for (Edge it : edges) {
                        int key_1 = DFS_CODE.get(aRmPath).from;
                        NavigableMap<Integer, NavigableMap<Integer, Projected>> root_1 = new_fwd_root.computeIfAbsent(key_1, k -> new TreeMap<>());
                        int key_2 = it.eLabel;
                        NavigableMap<Integer, Projected> root_2 = root_1.computeIfAbsent(key_2, k -> new TreeMap<>());
                        int key_3 = TRANS.get(id).get(it.to).label;
                        Projected root_3 = root_2.get(key_3);
                        if (root_3 == null) {
                            root_3 = new Projected();
                            root_2.put(key_3, root_3);
                        }
                        root_3.push(id, it, aProjected);
                    }
        }

        int Benefitscore_max_backward = -1;
        Integer maxToc_max_backward = 0;
        Integer index1_max_backward = 0;
        Integer eLabel_max_backward = 0;
        Projected projected_max_backward = null;
        
        
        // Test all extended substructures.
        // backward
        for (Entry<Integer, NavigableMap<Integer, Projected>> to : new_bck_root.entrySet()) {
            for (Entry<Integer, Projected> eLabel : to.getValue().entrySet()) {
            	
                int benefitscore =  getBenefitScore_Initial(eLabel.getValue());
                //System.out.println("benefitscore: " +benefitscore);
                if(benefitscore >= Benefitscore_max_backward) {
                	Benefitscore_max_backward = benefitscore;
                	maxToc_max_backward = maxToc;
                	index1_max_backward = to.getKey();
                	eLabel_max_backward = eLabel.getKey();
                	projected_max_backward = eLabel.getValue();
                }
            	//DFS_CODE.push(maxToc, to.getKey(), -1, eLabel.getKey(), -1);
            	//project_Initial(eLabel.getValue());
                //DFS_CODE.pop();
            }
        }
        
        
        int Benefitscore_max_forward = -1;
        Integer from_max_forward = 0;
        Integer maxToc_max_forward = 0;
        Integer index1_max_forward = 0;
        Integer eLabel_max_forward = 0;
        Projected projected_max_forward = null;
        //note:  There are many forward edges, so they should be visited descendingly resp.  from_vertex.
        //            e.g., let rightmost path be (1,2,4), then forward edges should be, (4,5), (2,5) and (1,5)
        // forward
        for (Entry<Integer, NavigableMap<Integer, NavigableMap<Integer, Projected>>> from : new_fwd_root.descendingMap().entrySet()) {
            for (Entry<Integer, NavigableMap<Integer, Projected>> eLabel : from.getValue().entrySet()) {
                for (Entry<Integer, Projected> toLabel : eLabel.getValue().entrySet()) {
                	
                    int benefitscore =  getBenefitScore_Initial(toLabel.getValue());
                    if(benefitscore>=Benefitscore_max_forward) {
                    	Benefitscore_max_forward = benefitscore;
                    	from_max_forward   = from.getKey();
                    	maxToc_max_forward = maxToc + 1;
                    	index1_max_forward = eLabel.getKey();
                    	eLabel_max_forward = toLabel.getKey();
                    	projected_max_forward = toLabel.getValue();
                    }
                	//DFS_CODE.push(from.getKey(), maxToc + 1, -1, eLabel.getKey(), toLabel.getKey());
                   // project(toLabel.getValue());
                   // DFS_CODE.pop();
                }
            }
        }
        //System.out.println("CurrentBenefitScore, Benefitscore_max_forward, Benefitscore_max_backward: "+ CurrentBenefitScore + "," + Benefitscore_max_forward + "," + Benefitscore_max_backward);
    	
        if(Benefitscore_max_forward >= CurrentBenefitScore && Benefitscore_max_forward>= Benefitscore_max_backward) {
        	//System.out.println("here1");
        	DFS_CODE.push(from_max_forward, maxToc_max_forward, -1, index1_max_forward, eLabel_max_forward);
        	project_Initial(projected_max_forward);
            DFS_CODE.pop();
        }else if(Benefitscore_max_backward >= CurrentBenefitScore &&  Benefitscore_max_backward >= Benefitscore_max_forward) {
        	//System.out.println("here2");
        	DFS_CODE.push(maxToc_max_backward, index1_max_backward, -1, eLabel_max_backward, -1);
        	project_Initial(projected_max_backward);
            DFS_CODE.pop();
        }else {
        	//System.out.println("here3");
        	//System.out.println("CurrentBenefitScore, Benefitscore_max_forward, Benefitscore_max_backward: "+ CurrentBenefitScore + "," + Benefitscore_max_forward + "," + Benefitscore_max_backward);
        	reportwithCoveredEdges_Initial(0, projected);
        	return ;
        }
    }

    private int support(Projected projected) {
        int oid = 0xffffffff;
        int size = 0;

        for (PDFS cur : projected) {
            if (oid != cur.id) {
                ++size;
            }
            oid = cur.id;
        }

        return size;
    }

    private boolean isMin() {
        if (DFS_CODE.size() == 1)
            return (true);

        DFS_CODE.toGraph(GRAPH_IS_MIN);
        DFS_CODE_IS_MIN.clear();
        
        // note: [vertex1.label][eLabel][vertex2.label] = Projected
        NavigableMap<Integer, NavigableMap<Integer, NavigableMap<Integer, Projected>>> root = new TreeMap<>();
        ArrayList<Edge> edges = new ArrayList<>();

        for (int from = 0; from < GRAPH_IS_MIN.size(); ++from)
            if (Misc.getForwardRoot(GRAPH_IS_MIN, GRAPH_IS_MIN.get(from), edges))
                for (Edge it : edges) {
                    int key_1 = GRAPH_IS_MIN.get(from).label;
                    NavigableMap<Integer, NavigableMap<Integer, Projected>> root_1 = root.computeIfAbsent(key_1, k -> new TreeMap<>());
                    int key_2 = it.eLabel;
                    NavigableMap<Integer, Projected> root_2 = root_1.computeIfAbsent(key_2, k -> new TreeMap<>());
                    int key_3 = GRAPH_IS_MIN.get(it.to).label;
                    Projected root_3 = root_2.get(key_3);
                    if (root_3 == null) {
                        root_3 = new Projected();
                        root_2.put(key_3, root_3);
                    }
                    // note: [vertex1.label][eLabel][vertex2.label] = Projected, but here the graph id is fixed as 0, since Projected is only for GRAPH_IS_MIN
                    root_3.push(0, it, null);
                }

        Entry<Integer, NavigableMap<Integer, NavigableMap<Integer, Projected>>> fromLabel = root.firstEntry();
        Entry<Integer, NavigableMap<Integer, Projected>> eLabel = fromLabel.getValue().firstEntry();
        Entry<Integer, Projected> toLabel = eLabel.getValue().firstEntry();
        // note: select the minimum edge as the as first DFS of DFS_CODE_IS_MIN
        DFS_CODE_IS_MIN.push(0, 1, fromLabel.getKey(), eLabel.getKey(), toLabel.getKey());

        return isMinProject(toLabel.getValue());
    }

    ///// note: similar to the project function
    private boolean isMinProject(Projected projected) {
    	//note: rmPath contains the edges (DFS) in the rightmost path. [i] = i-th DFS.  It is generated from the last DFS to the first DFS 
        ArrayList<Integer> rmPath = DFS_CODE_IS_MIN.buildRMPath(); 
        
        //note: fromlabel of first DFS in  DFS_CODE_IS_MIN
        int minLabel = DFS_CODE_IS_MIN.get(0).fromLabel;
        
        //note: rightmost vertex in  DFS_CODE_IS_MIN
        int maxToc = DFS_CODE_IS_MIN.get(rmPath.get(0)).to;

        {
            NavigableMap<Integer, Projected> root = new TreeMap<>();
            boolean flg = false;
            int newTo = 0;

            for (int i = rmPath.size() - 1; !flg && i >= 1; --i) {
                for (PDFS cur : projected) {
                    History history = new History(GRAPH_IS_MIN, cur);
                    Edge e = Misc.getBackward(GRAPH_IS_MIN, history.get(rmPath.get(i)), history.get(rmPath.get(0)),
                            history);
                    if (e != null) {
                        int key_1 = e.eLabel;
                        Projected root_1 = root.get(key_1);
                        if (root_1 == null) {
                            root_1 = new Projected();
                            root.put(key_1, root_1);
                        }
                        root_1.push(0, e, cur);
                        newTo = DFS_CODE_IS_MIN.get(rmPath.get(i)).from;
                        flg = true;
                    }
                }
            }

            if (flg) {
                Entry<Integer, Projected> eLabel = root.firstEntry();
                DFS_CODE_IS_MIN.push(maxToc, newTo, -1, eLabel.getKey(), -1);
                if (DFS_CODE.get(DFS_CODE_IS_MIN.size() - 1)
                        .notEqual(DFS_CODE_IS_MIN.get(DFS_CODE_IS_MIN.size() - 1)))
                    return false;
                return isMinProject(eLabel.getValue());
            }
        }

        {
            boolean flg = false;
            int newFrom = 0;
            NavigableMap<Integer, NavigableMap<Integer, Projected>> root = new TreeMap<>();
            ArrayList<Edge> edges = new ArrayList<>();

            for (PDFS cur : projected) {
                History history = new History(GRAPH_IS_MIN, cur);
                if (Misc.getForwardPure(GRAPH_IS_MIN, history.get(rmPath.get(0)), minLabel, history, edges)) {
                    flg = true;
                    newFrom = maxToc;
                    for (Edge it : edges) {
                        int key_1 = it.eLabel;
                        NavigableMap<Integer, Projected> root_1 = root.computeIfAbsent(key_1, k -> new TreeMap<>());
                        int key_2 = GRAPH_IS_MIN.get(it.to).label;
                        Projected root_2 = root_1.get(key_2);
                        if (root_2 == null) {
                            root_2 = new Projected();
                            root_1.put(key_2, root_2);
                        }
                        root_2.push(0, it, cur);
                    }
                }
            }

            for (int i = 0; !flg && i < rmPath.size(); ++i) {
                for (PDFS cur : projected) {
                    History history = new History(GRAPH_IS_MIN, cur);
                    if (Misc.getForwardRmPath(GRAPH_IS_MIN, history.get(rmPath.get(i)), minLabel, history, edges)) {
                        flg = true;
                        newFrom = DFS_CODE_IS_MIN.get(rmPath.get(i)).from;
                        for (Edge it : edges) {
                            int key_1 = it.eLabel;
                            NavigableMap<Integer, Projected> root_1 = root.computeIfAbsent(key_1, k -> new TreeMap<>());
                            int key_2 = GRAPH_IS_MIN.get(it.to).label;
                            Projected root_2 = root_1.get(key_2);
                            if (root_2 == null) {
                                root_2 = new Projected();
                                root_1.put(key_2, root_2);
                            }
                            root_2.push(0, it, cur);
                        }
                    }
                }
            }

            if (flg) {
                Entry<Integer, NavigableMap<Integer, Projected>> eLabel = root.firstEntry();
                Entry<Integer, Projected> toLabel = eLabel.getValue().firstEntry();
                DFS_CODE_IS_MIN.push(newFrom, maxToc + 1, -1, eLabel.getKey(), toLabel.getKey());
                if (DFS_CODE.get(DFS_CODE_IS_MIN.size() - 1)
                        .notEqual(DFS_CODE_IS_MIN.get(DFS_CODE_IS_MIN.size() - 1)))
                    return false;
                return isMinProject(toLabel.getValue());
            }
        }

        return true;
    }
}
