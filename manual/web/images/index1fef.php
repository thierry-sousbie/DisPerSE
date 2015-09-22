<?xml version="1.0" encoding="utf-8"?><feed xmlns="http://www.w3.org/2005/Atom"
  xmlns:dc="http://purl.org/dc/elements/1.1/"
  xmlns:wfw="http://wellformedweb.org/CommentAPI/"
  xml:lang="en">
  
  <title type="html">DisPerSE - persistent structures identification - Tutorial</title>
  <subtitle type="html">Automatic identification of persistent structures in 2D or 3D.
keywords: Morse complex, topology, peak, void, source, wall, filament, cosmic web, cosmology, structure identification.</subtitle>
  <link href="http://localhost/dotclear/index.php?feed/category/Quick-start/atom" rel="self" type="application/atom+xml"/>
  <link href="http://localhost/dotclear/index.php?" rel="alternate" type="text/html"
  title="Automatic identification of persistent structures in 2D or 3D.
keywords: Morse complex, topology, peak, void, source, wall, filament, cosmic web, cosmology, structure identification."/>
  <updated>2013-01-24T03:58:22+01:00</updated>
  <author>
    <name>thierry sousbie</name>
  </author>
  <id>urn:md5:f62b3bf39b9948523836119476648bbf</id>
  <generator uri="http://www.dotclear.org/">Dotclear</generator>
  
    
  <entry>
    <title>Point sample: 3D walls and filaments</title>
    <link href="http://localhost/dotclear/index.php?post/Example-2" rel="alternate" type="text/html"
    title="Point sample: 3D walls and filaments" />
    <id>urn:md5:94d68b84101013d641d848c9a2d563b2</id>
    <published>2320-05-09T06:59:00+01:00</published>
    <updated>2012-10-02T17:01:36+02:00</updated>
    <author><name>thierry sousbie</name></author>
        <dc:subject>Tutorial</dc:subject>
            
    <content type="html">    &lt;p&gt;This tutorial follows the 2D discrete case example introduced in the &lt;a href=&quot;http://localhost/dotclear/index.php?post/Example-1&quot;&gt;previous section&lt;/a&gt; so you should probably read it first if you haven't done so yet.
&lt;br /&gt;&lt;br /&gt;
We will compute the filaments and walls from the downsampled dark matter N-body simulation snapshot &lt;em&gt;${DISPERSE}/data/simu_32_id.gad&lt;/em&gt;. The distribution has periodic boundary conditions and we therefore start by computing its Delaunay tessellation with the following command:
&lt;br /&gt;&lt;br /&gt;&lt;/p&gt;
&lt;div style=&quot;text-align: center; border: 1px dotted gray;&quot;&gt;delaunay_3D simu_32_id.gad -periodic&lt;/div&gt;


&lt;p&gt;&lt;br /&gt;
The program should output a file called &lt;em&gt;simu_32_id.gad.NDnet&lt;/em&gt; that can be used directly with &lt;a href=&quot;http://localhost/dotclear/index.php?post/mse&quot;&gt;mse&lt;/a&gt; to compute the filaments and walls of the distribution. In 3D, walls are represented by the &lt;a href=&quot;http://localhost/dotclear/index.php?post/definitions&quot;&gt;ascending 2-manifolds&lt;/a&gt; that originate from critical points of critical index &lt;em&gt;1&lt;/em&gt; so they can be computed with option &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#dumpmanifolds&quot;&gt;-dumpManifolds&lt;/a&gt; JE1a&lt;/em&gt; (letter &lt;em&gt;E&lt;/em&gt; stands for &lt;em&gt;Extended&lt;/em&gt;, which means that the ascending 1-manifolds and ascending 0-manifolds  at the boundary of the walls will also be added, see option &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#dumpmanifolds&quot;&gt;-dumpManifolds&lt;/a&gt;&lt;/em&gt; of &lt;a href=&quot;http://localhost/dotclear/index.php?post/mse&quot;&gt;mse&lt;/a&gt;) .Note that by default, subsets of the walls that are common to several walls are merged together, so you may want to try &lt;em&gt;-dumpManifolds JP1a&lt;/em&gt; instead if you are interested in studying the properties of individual walls. The filaments are computed with option &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#upskl&quot;&gt;-upSkl&lt;/a&gt;&lt;/em&gt; or &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#dumparcs&quot;&gt;-dumparcs&lt;/a&gt; U&lt;/em&gt; which are equivalent. Using interactive mode to select a persistence threshold (use &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#nsig&quot;&gt;-nsig&lt;/a&gt; 3.5&lt;/em&gt; instead of &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#interactive&quot;&gt;-interactive&lt;/a&gt;&lt;/em&gt; if &lt;a href=&quot;http://localhost/dotclear/index.php?post/pdview&quot;&gt;pdview&lt;/a&gt; could not compile on your computer), we run the following command:
&lt;br /&gt;&lt;br /&gt;&lt;/p&gt;
&lt;div style=&quot;text-align: center; border: 1px dotted gray;&quot;&gt;mse simu_32_id.gad.NDnet -dumpManifolds JE1a -upSkl -interactive&lt;/div&gt;


&lt;p&gt;&lt;br /&gt;
which after some quick computations, should display the following persistence diagram.
&lt;br /&gt;&lt;br /&gt;
&lt;img src=&quot;http://localhost/dotclear/public/per_diag_3D.jpg&quot; alt=&quot;per_diag_3D.jpg&quot; style=&quot;display:block; margin:0 auto;&quot; /&gt;
&lt;br /&gt;
On this diagram, each type of pair is represented by a mark of different color : blue for minima/1-saddle, green for 1-saddle/2-saddle and red for 2-saddle/maxima. Note that the blue dots in the lower left corner stand for the voids in the distribution, so it is expected that we find only few of them compared to the red dots which represent collapsed structure (the simulation is &lt;em&gt;50 Mpc&lt;/em&gt; large). Each point on the diagram representing a topological feature, we would like to select only those points that stand out of the general distribution in terms of persistence in order to keep only the most meaningful structures (the persistence selection threshold is set by clicking while holding &lt;em&gt;Ctrl&lt;/em&gt; key). The diagram suggests that persistence threshold should be above &lt;em&gt;~3-sigmas&lt;/em&gt; so we select a threshold of &lt;em&gt;~3.5 sigmas&lt;/em&gt; and click on the &lt;em&gt;done&lt;/em&gt; button to confirm the selection. The program continues and should output the filaments in file &lt;em&gt;simu_32_id.gad.NDnet_s3.52.up.NDskl&lt;/em&gt; and the walls in file &lt;em&gt;simu_32_id.gad.NDnet_s3.52_manifolds_JE1a.NDnet&lt;/em&gt;.
&lt;br /&gt;&lt;br /&gt;
It is instructive at that point to read &lt;a href=&quot;http://localhost/dotclear/index.php?post/mse&quot;&gt;mse&lt;/a&gt; output, and in particular the following section :
&lt;br /&gt;&lt;br /&gt;&lt;/p&gt;

&lt;pre&gt;****** Simplifying complex ******
Starting Morse-Smale complex simplification.
Computing persistence pairs ... SKIPPED.
Sampling noise level was set to 3.5-sigma.
Cancelling pairs with persistance ratio &amp;lt; [2.66e+00,3.14e+00,4.42e+01].
   Tagging arcs ... done.
   Cancelling pairs (smart) ... (4140 rem.)
done.
   Cancellation took 0.10s (4140 canceled, 1 conflicts, 10 loops).&lt;/pre&gt;

&lt;p&gt;&lt;br /&gt;&lt;/p&gt;


&lt;p&gt;what &lt;em&gt;mse&lt;/em&gt; tells us here is that it successfully removed &lt;em&gt;4140&lt;/em&gt; persistence pairs, but was not able to cancel &lt;em&gt;11&lt;/em&gt; more although their persistence was lower than the threshold. The &lt;em&gt;10 loops&lt;/em&gt; represent persistence pairs that were connected by more than one arc at the moment of their cancellation. Topologically speaking, their cancellation is not a valid operation as it would result in a discrete gradient loop (integral lines cannot form loops in &lt;a href=&quot;http://localhost/dotclear/index.php?post/general-concepts&quot;&gt;Morse-theory&lt;/a&gt;). Note that the existence of non-cancellable pairs is not a bug of DisPerSE but rather a feature of Morse theory. However, not cancelling them may result in a few spurious non-significant pieces of filaments remaining in the ouput file (this is very problematic though as the impact on the identified filaments is very small in general) so one can use option &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#forceloops&quot;&gt;-forceLoops&lt;/a&gt;&lt;/em&gt; to remove them anyway (which is perfectly fine if your goal is only to identify structures). Using this option will also solve the &lt;em&gt;1 conflict&lt;/em&gt; identified by mse, as conflicts are the result of a non-cancellable pair blocking the cancellation of a valid cancellable pair.
&lt;br /&gt;&lt;br /&gt;
We can now rerun &lt;em&gt;mse&lt;/em&gt; while skipping the computations by loading the backup file &lt;em&gt;simu_32_id.gad.NDnet.MSC&lt;/em&gt; that mse generated on the previous run with option &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#loadmsc&quot;&gt;-loadMSC&lt;/a&gt;&lt;/em&gt;:
&lt;br /&gt;&lt;br /&gt;&lt;/p&gt;
&lt;div style=&quot;text-align: center; border: 1px dotted gray;&quot;&gt;mse simu_32_id.gad.NDnet -dumpManifolds JE1a -upSkl -forceLoops -loadMSC simu_32_id.gad.NDnet.MSC  -nsig 3.52&lt;/div&gt;


&lt;p&gt;&lt;br /&gt;
and the ouput should contain the following line:
&lt;br /&gt;&lt;br /&gt;&lt;/p&gt;

&lt;pre&gt;   Cancellation took 0.10s (4151 canceled, 0 conflicts, 11 forced loops).&lt;/pre&gt;

&lt;p&gt;&lt;br /&gt;
indicating that all the &lt;em&gt;4151&lt;/em&gt; non persistent pairs could be cancelled.
&lt;br /&gt;&lt;br /&gt;
We can now smooth and convert the resulting files to a more portable format:
&lt;br /&gt;&lt;br /&gt;&lt;/p&gt;

&lt;pre&gt;skelconv simu_32_id.gad.NDnet_s3.52.up.NDskl -smooth 10 -to vtp
netconv simu_32_id.gad.NDnet -to vtu
netconv simu_32_id.gad.NDnet_s3.52_manifolds_JE1a.NDnet -smooth 10 -to vtu&lt;/pre&gt;

&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;walls&quot;&gt;&lt;/a&gt;
and visualize the result with &lt;a href=&quot;http://www.paraview.org/&quot; title=&quot;paraview&quot;&gt;paraview&lt;/a&gt;. The configuration file &lt;em&gt;${DISPERSE}/data/simu_32_per.pvsm&lt;/em&gt; can be directly loaded in paraview (&lt;em&gt;File-&amp;gt;Load State&lt;/em&gt;) to make the following plot: (note that it is a bit tricky to deal with periodic boundary conditions when visualizing, clipping the boundaries is usually the easiest way)
&lt;br /&gt;&lt;br /&gt;
&lt;img src=&quot;http://localhost/dotclear/public/simu_3D_result_small.jpg&quot; alt=&quot;simu_3D_result_small.jpg&quot; style=&quot;display:block; margin:0 auto;&quot; /&gt;
&lt;br /&gt;&lt;br /&gt;
&lt;strong&gt;Note&lt;/strong&gt;: contrary to the previous example where we computed 2D voids (see &lt;a href=&quot;http://localhost/dotclear/index.php?post/Example-1#voids&quot;&gt;here&lt;/a&gt;), it is problematic in this case to identify each individual wall and paint it with a specific color. Indeed, contrary to voids, the 2-manifolds that represent walls can overlap, and so a unique &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/networks#source&quot;&gt;source_index&lt;/a&gt;&lt;/em&gt; representing the critical point from which each wall originates cannot be assigned to each simplex (i.e. in that case, triangle). A solution to this problem would have been to specify that we wanted identical simplices to be stored as many times as they appear in different walls even though they are the same. This could have been achieved with flag &lt;em&gt;P&lt;/em&gt; of option &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#dumpmanifolds&quot;&gt;-dumpManifolds&lt;/a&gt;&lt;/em&gt; in &lt;em&gt;mse&lt;/em&gt;:
&lt;br /&gt;&lt;br /&gt;&lt;/p&gt;
&lt;div style=&quot;text-align: center; border: 1px dotted gray;&quot;&gt;mse simu_32_id.gad.NDnet -dumpManifolds JEP1a -upSkl -forceLoops -loadMSC simu_32_id.gad.NDnet.MSC  -nsig 3.52&lt;/div&gt;


&lt;p&gt;&lt;br /&gt;
As a result, the output file is larger, but it would also allow cross-linking the information contained in the geometry of the walls and that in the skeleton files or persistence pairs network (see this &lt;a href=&quot;http://localhost/dotclear/index.php?post/Example-1#voids&quot;&gt;note&lt;/a&gt;).&lt;/p&gt;</content>
    
    

    
      </entry>
    
  <entry>
    <title>Point sample: 2D voids and filaments</title>
    <link href="http://localhost/dotclear/index.php?post/Example-1" rel="alternate" type="text/html"
    title="Point sample: 2D voids and filaments" />
    <id>urn:md5:78ae47fc599db688f2e112bfcca6dce2</id>
    <published>2310-05-09T06:59:00+01:00</published>
    <updated>2012-10-02T16:55:23+02:00</updated>
    <author><name>thierry sousbie</name></author>
        <dc:subject>Tutorial</dc:subject>
            
    <content type="html">    &lt;p&gt;This tutorial uses file  &lt;em&gt;${DISPERSE}/data/simu_2D.ND&lt;/em&gt; which contains the projected distribution of the particles in a slice of an N-Body dark matter simulation of a &lt;em&gt;50 Mpc&lt;/em&gt; large chunk of the universe. Note that although we use a 2D slice for convenience here, the following procedure would also work with the full 3D distribution.
&lt;br /&gt;
&lt;br /&gt;
In this example, the voids and filaments of the particle distribution will be computed from the underlying density function traced by the particles, so we first need to define the topology of the space (i.e. we need a notion of neighborhood for each particle) and to estimate the density function. This is achieved using the Delaunay tessellation of the particle distribution, which can be computed with &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/Usage&quot;&gt;delaunay_2D&lt;/a&gt;&lt;/em&gt;.
&lt;br /&gt;
&lt;br /&gt;
The distribution of the particles is periodic, so we could simply run:
&lt;br /&gt;&lt;br /&gt;&lt;/p&gt;
&lt;div style=&quot;text-align: center; border: 1px dotted gray;&quot;&gt;delaunay_2D simu_2D.ND -periodic&lt;/div&gt;


&lt;p&gt;&lt;br /&gt;
However, in this example we are going to pretend that the distribution has boundaries, which will make visualization easier (if boundary conditions are periodic, a filament crossing a boundary reappears on the other side, causing visual artifacts). We nevertheless would like to use the fact that boundary conditions are indeed periodic to correctly estimate density on the boundaries, which can be achieved using option &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/Usage#btype&quot;&gt;-btype&lt;/a&gt;&lt;/em&gt;:
&lt;br /&gt;&lt;br /&gt;&lt;/p&gt;
&lt;div style=&quot;text-align: center; border: 1px dotted gray;&quot;&gt;delaunay_2D simu_2D.ND -btype periodic&lt;/div&gt;


&lt;p&gt;&lt;br /&gt;
The output file is named &lt;em&gt;simu_2D.ND.NDnet&lt;/em&gt; by default, and it contains the Delauany tessellation as an &lt;a href=&quot;http://localhost/dotclear/index.php?post/networks&quot;&gt;unstructured network&lt;/a&gt; in &lt;a href=&quot;http://localhost/dotclear/index.php?post/NDnet-format&quot;&gt;NDnet&lt;/a&gt; format, as well as the density function estimated using DTFE at each vertex (i.e. the additional vertex data labeled &lt;em&gt;field_value&lt;/em&gt; in the NDnet file). This file can be fed directly to &lt;a href=&quot;http://localhost/dotclear/index.php?post/mse&quot;&gt;mse&lt;/a&gt;, but we have to choose a persistence threshold first.
&lt;br /&gt;
&lt;br /&gt;
Because the input file is a Delaunay tessellation, persistence is expressed as a ratio (density has a close to lognormal PDF) and the persistence threshold as a &lt;em&gt;number of sigmas&lt;/em&gt; representing the probability that a persistence pair with given persistence ratio happens in a pure random discrete Poisson distribution. It is usually safe in such case to select a threshold between 2-sigmas and 4-sigmas, and we could directly compute the filaments of the distribution with the following command (see &lt;a href=&quot;http://localhost/dotclear/index.php?post/mse&quot;&gt;mse&lt;/a&gt; for more info):
&lt;br /&gt;&lt;br /&gt;&lt;/p&gt;
&lt;div style=&quot;text-align: center; border: 1px dotted gray;&quot;&gt;mse simu_2D.ND.NDnet -nsig 3 -upSkl&lt;/div&gt;


&lt;p&gt;&lt;br /&gt;
However, in this example, we will use a persistence diagram to decide the correct threshold (in case you do not have &lt;a href=&quot;http://localhost/dotclear/index.php?post/pdview&quot;&gt;pdview&lt;/a&gt;, replace &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#interactive&quot;&gt;-interactive&lt;/a&gt;&lt;/em&gt; with &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#nsig&quot;&gt;-nsig 3&lt;/a&gt;&lt;/em&gt; in the following commannd) and we would also like to identify the voids in the distribution (i.e. the ascending 3-manifolds originating from critical points of critical index 0; see option &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#dumpmanifolds&quot;&gt;-dumpManifolds&lt;/a&gt;&lt;/em&gt; of &lt;a href=&quot;http://localhost/dotclear/index.php?post/mse&quot;&gt;mse&lt;/a&gt;):
&lt;br /&gt;&lt;br /&gt;&lt;/p&gt;
&lt;div style=&quot;text-align: center; border: 1px dotted gray;&quot;&gt;mse simu_2D.ND.NDnet -interactive -upSkl -dumpManifolds J0a&lt;/div&gt;


&lt;p&gt;&lt;br /&gt;
After running this command line, you should be able to visualize the following diagrams:
&lt;br /&gt;&lt;br /&gt;
&lt;img src=&quot;http://localhost/dotclear/public/per_diag_small.jpg&quot; alt=&quot;per_diag_small.jpg&quot; style=&quot;display:block; margin:0 auto;&quot; /&gt;
&lt;br /&gt;
On this picture, the maxima/saddle persistence pairs diagram is represented on the left and the minima/saddle one on the right. The fact that the second type of pairs are clearly less persistent reflects the known fact that dark matter haloes are very prominent structures because they belong to regions in the non linear regime where density increases exponentially fast. One nice feature of persistence diagrams is that they make it easy to choose the threshold above which underlying topological structures stand out of the noise : in the present case, we want to select the high persistence tails only that clearly stand out from the bunch of pairs generated by noise. This gives us a threshold of around &lt;em&gt;3-sigmas&lt;/em&gt; that can be visually selected by clicking while holding &lt;em&gt;Ctrl&lt;/em&gt; key down. Once the threshold is selected, simply press &lt;em&gt;Done&lt;/em&gt; button.
&lt;br /&gt;
&lt;br /&gt;
After this operation, &lt;a href=&quot;http://localhost/dotclear/index.php?post/mse&quot;&gt;mse&lt;/a&gt; should output two files: the filaments in a &lt;a href=&quot;http://localhost/dotclear/index.php?post/skeleton-formats&quot;&gt;skeleton type&lt;/a&gt; file (&lt;em&gt;simu_2D.ND.NDnet_s3.up.NDskl&lt;/em&gt;) and the voids in a &lt;a href=&quot;http://localhost/dotclear/index.php?post/networks&quot;&gt;network type&lt;/a&gt; file (&lt;em&gt;simu_2D.ND.NDnet_s3_manifolds_J0a.NDnet&lt;/em&gt;).
&lt;br /&gt;
&lt;br /&gt;
By default, &lt;em&gt;mse&lt;/em&gt; associates maxima to vertices while minima are associated to d-&lt;a href=&quot;http://localhost/dotclear/index.php?post/definitions&quot;&gt;simplices&lt;/a&gt; (where d is the dimension of space). As a result, the network representing the voids defines a set of triangles each associated to a given void. This may not be very practical for post-treatment or visualization, but we can use option &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#vertexasminima&quot;&gt;-vertexAsMinima&lt;/a&gt;&lt;/em&gt; to reverse the convention (this will overwrite the previous file if the threshold is identical):
&lt;br /&gt;&lt;br /&gt;&lt;/p&gt;
&lt;div style=&quot;text-align: center; border: 1px dotted gray;&quot;&gt;mse simu_2D.ND.NDnet -nsig 3 -vertexAsMinima -dumpManifolds J0a&lt;/div&gt;


&lt;p&gt;&lt;br /&gt;
We can now visualize the result by converting the files to vtk formats using &lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv&quot;&gt;skelconv&lt;/a&gt; and &lt;a href=&quot;http://localhost/dotclear/index.php?post/netconv&quot;&gt;netconv&lt;/a&gt;. The skeleton is smoothed and converted to &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/vtk-formats&quot;&gt;vtp&lt;/a&gt;&lt;/em&gt; format:
&lt;br /&gt;&lt;br /&gt;&lt;/p&gt;
&lt;div style=&quot;text-align: center; border: 1px dotted gray;&quot;&gt;skelconv simu_2D.ND.NDnet_s3.up.NDskl -smooth 10 -to vtp&lt;/div&gt;


&lt;p&gt;&lt;br /&gt;
and the input tesselation and output voids are converted to &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/vtk-network-format&quot;&gt;vtu&lt;/a&gt;&lt;/em&gt; format:
&lt;br /&gt;&lt;br /&gt;&lt;/p&gt;
&lt;div style=&quot;text-align: center; border: 1px dotted gray;&quot;&gt;
netconv simu_2D.ND.NDnet_s3_manifolds_J0a.NDnet -to vtu&lt;/br&gt;
netconv simu_2D.ND.NDnet -to vtu
&lt;/div&gt;


&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;voids&quot;&gt;&lt;/a&gt;
the result should look as follows when visualized with &lt;em&gt;&lt;a href=&quot;http://www.paraview.org/&quot; title=&quot;paraview&quot;&gt;paraview&lt;/a&gt;&lt;/em&gt;: (the configuration file &lt;em&gt;${DISPERSE}/data/simu_2D_nonper.pvsm&lt;/em&gt; can be directly loaded in paraview (&lt;em&gt;File-&amp;gt;Load State&lt;/em&gt;) to produce the following plot)
&lt;br /&gt;&lt;br /&gt;
&lt;img src=&quot;http://localhost/dotclear/public/simu_2D_tuto_small.jpg&quot; alt=&quot;simu_2D_tuto_small.jpg&quot; style=&quot;display:block; margin:0 auto;&quot; /&gt;
&lt;br /&gt;&lt;br /&gt;
&lt;strong&gt;Note&lt;/strong&gt;: on the right figure, a different color is associated to the &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/networks#source&quot;&gt;source_index&lt;/a&gt;&lt;/em&gt; of each particle depending on which void it belongs to: each void originates from a different minimum of the field and each critical point has a specific index. This &lt;em&gt;source_index&lt;/em&gt; corresponds to the  position of this critical point as a node in skeleton files or as a vertex in persistence pair networks obtained with &lt;a href=&quot;http://localhost/dotclear/index.php?post/mse&quot;&gt;mse&lt;/a&gt; (options &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#dumparcs&quot;&gt;-dumpArcs&lt;/a&gt;&lt;/em&gt; and &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#ppairs&quot;&gt;-ppairs&lt;/a&gt;&lt;/em&gt;). It can therefore also be used to cross-link information contained in manifold networks and skeletons or persistence pairs, such as for instance retrieving the hierarchy of the voids (i.e. using the &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/networks#parent&quot;&gt;parent_index&lt;/a&gt;&lt;/em&gt; of the critical points in &lt;a href=&quot;http://localhost/dotclear/index.php?post/skeleton-formats#parent&quot;&gt;skeleton files&lt;/a&gt; and persistence pairs networks).&lt;/p&gt;</content>
    
    

    
      </entry>
    
  <entry>
    <title>Persistence diagrams and filaments</title>
    <link href="http://localhost/dotclear/index.php?post/regular-grid-%3A-filaments" rel="alternate" type="text/html"
    title="Persistence diagrams and filaments" />
    <id>urn:md5:e4ef47e75cc7c98d412cb709fd81079a</id>
    <published>2305-05-09T06:59:00+01:00</published>
    <updated>2012-07-11T20:36:28+02:00</updated>
    <author><name>thierry sousbie</name></author>
        <dc:subject>Tutorial</dc:subject>
            
    <content type="html">    &lt;p&gt;In this tutorial, we will identify the persistent filamentary structures in a diHydrogen column density map of the interstellar medium in the Eagle Nebula (M16) as observed by &lt;a href=&quot;http://www.esa.int/herschel&quot; title=&quot;herschel&quot;&gt;Herschel&lt;/a&gt; . The map is encoded in a 2D double precision FITS image, with pixels set to &lt;em&gt;NaN&lt;/em&gt; in unobserved regions:
&lt;br /&gt;&lt;br /&gt;
&lt;img src=&quot;http://localhost/dotclear/public/M16_blue.png&quot; alt=&quot;M16_blue.png&quot; style=&quot;display:block; margin:0 auto;&quot; /&gt;
&lt;br /&gt;&lt;br /&gt;
By default, &lt;a href=&quot;http://localhost/dotclear/index.php?post/mse&quot;&gt;mse&lt;/a&gt; will automatically mask NaN valued pixels so the FITS file could theoretically be fed directly to it. However, because of the complexity of the geometry of the observed region, it is preferable to use a mask to limit the analysis to a well defined region of interest (in particular, it is preferable that this region is simply connected - i.e. does not contain holes -). Designing a mask is relatively simple, as it consists in an image with pixels set to 0 in non-masked region and any other value in the masked regions (this behavior can be changed by adding a trailing '~' to the mask filename, see option &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#mask&quot;&gt;-mask&lt;/a&gt;&lt;/em&gt;). The following mask was created with &lt;em&gt;&lt;a href=&quot;http://www.gimp.org/&quot; title=&quot;gimp&quot;&gt;Gimp&lt;/a&gt;&lt;/em&gt; and saved as a &lt;em&gt;.PNG&lt;/em&gt; file (any readable &lt;a href=&quot;http://localhost/dotclear/index.php?post/field-formats&quot;&gt;field file format&lt;/a&gt; is acceptable):
&lt;br /&gt;&lt;br /&gt;
&lt;img src=&quot;http://localhost/dotclear/public/M16_mask.png&quot; alt=&quot;M16_mask.png&quot; style=&quot;display:block; margin:0 auto;&quot; /&gt;
&lt;br /&gt;&lt;br /&gt;
&lt;br /&gt;&lt;br /&gt;
With this mask, we can now use &lt;em&gt;mse&lt;/em&gt; to identify the filamentary structures in the map (option &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#upskl&quot;&gt;-upSkl&lt;/a&gt;&lt;/em&gt;) with the help of an interactive persitence diagram (option &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#interactive&quot;&gt;-interactive&lt;/a&gt;&lt;/em&gt;) to select the correct persistent threshold as we do not want to impose any à-priori on its value. We therefore run the following command:
&lt;br /&gt;&lt;br /&gt;&lt;/p&gt;
&lt;div style=&quot;text-align: center; border: 1px dotted gray;&quot;&gt;mse m16_zo_no70um_density.fits -mask M16_mask.png -interactive -upSkl&lt;/div&gt;


&lt;p&gt;&lt;br /&gt;&lt;br /&gt;
After some computations, the following diagram should appear (use buttons &lt;em&gt;logX&lt;/em&gt;/&lt;em&gt;logY&lt;/em&gt; to switch logarithmic coordinate and buttons &lt;em&gt;1&lt;/em&gt;, &lt;em&gt;2&lt;/em&gt; and &lt;em&gt;3&lt;/em&gt; to show/hide different types of persistence pairs):
&lt;br /&gt;&lt;br /&gt;
&lt;img src=&quot;http://localhost/dotclear/public/M16_persistence_diagram_2.png&quot; alt=&quot;M16_persistence_diagram_2.png&quot; style=&quot;display:block; margin:0 auto;&quot; /&gt;
&lt;br /&gt;&lt;br /&gt;
On this diagram, each dot represents a &lt;a href=&quot;http://localhost/dotclear/index.php?post/definitions&quot;&gt;persistence pair&lt;/a&gt; of critical points (see &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/Persistence-and-simplification&quot;&gt;overview&lt;/a&gt;&lt;/em&gt; section), with its X coordinate being the value at the lowest critical point of the two (i.e. the &lt;em&gt;background&lt;/em&gt; density) and its Y coordinate the persitence of the pair (i.e. the value difference between the two, interpreted as &lt;em&gt;contrast&lt;/em&gt; of the topological feature it represents with respect to its background). Points located in the right side of this diagram represent features located in high valued regions of the map, while points in the upper side represent features that clearly stand out from their background. More specifically, the green squares stand for pairs of type &lt;em&gt;1&lt;/em&gt; (i.e. saddle-point / maximum pairs) while blue disks stand for pairs of type &lt;em&gt;0&lt;/em&gt; (i.e. minimum saddle-point pairs). This can be easily visualized on the following image, which represents a zoom on the map with critical points identified as blue squares, green crosses and red triangles for minima, saddle points and maxima respectively.  On this plot, the identified filaments are represented in black, and each persistence pair of type &lt;em&gt;0&lt;/em&gt; or &lt;em&gt;1&lt;/em&gt;  is represented as a blue or green segment linking its two critical points together (each segment corresponds to one point in the persistence diagram above):
&lt;br /&gt;&lt;br /&gt;
&lt;img src=&quot;http://localhost/dotclear/public/M16_skel_nocut_ppairs_crit.png&quot; alt=&quot;M16_skel_nocut_ppairs_crit.png&quot; style=&quot;display:block; margin:0 auto;&quot; /&gt;
&lt;br /&gt;&lt;br /&gt;&lt;/p&gt;


&lt;p&gt;Basically, in 2D, removing a pair of type &lt;em&gt;0&lt;/em&gt; (i.e. cancelling two critical points linked by a blue segment) amounts to deleting a piece of filament (i.e. 2 arcs departing from the saddle point) by merging two voids while removing a pair of type &lt;em&gt;1&lt;/em&gt; (i.e. cancelling two critical points linked by a green segment) amounts to gluing pieces of filaments into longer ones  (i.e. two arcs departing from the saddle point with one other arc departing from the maximum). When selecting the persistence threshold without à-priori, we therefore want to keep only the pairs that stand out of the general distribution on the Y-axis (i.e. with high persitence) so that only the most significant topological features remain. Selecting a threshold of &lt;em&gt;~4.2E21&lt;/em&gt; as shown on the persistence diagram above seems quite reasonable (using &lt;em&gt;Ctrl+left Click&lt;/em&gt;, and clicking on the &lt;em&gt;DONE&lt;/em&gt; button to proceed).
&lt;br /&gt;&lt;br /&gt;
&lt;a name=&quot;trim_exp&quot;&gt;&lt;/a&gt;
Visualizing a persistence diagram is a good way to determine the persistence threshold when assuming no previous knowledge of the data set. Note however that persistence could also have been set à-priori to the estimated amplitude of what is considered non-meaningful features or noise with option &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#cut&quot;&gt;-cut &amp;lt;val&amp;gt;&lt;/a&gt;&lt;/em&gt; instead of &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#interactive&quot;&gt;-interactive&lt;/a&gt;&lt;/em&gt;. It is nonetheless always useful to check the persistence diagram as it deeply reflects the topological nature of the distribution and can be used to understand it better. For instance, the following diagram (represented by a 2D histogram of the pairs distribution, instead of the pairs themselves) was computed from a totally different map, and it is a good example:
&lt;br /&gt;&lt;br /&gt;
&lt;img src=&quot;http://localhost/dotclear/public/ic5146_col_dens_pers_diag.png&quot; alt=&quot;ic5146_col_dens_pers_diag.png&quot; style=&quot;display:block; margin:0 auto;&quot; /&gt;
&lt;br /&gt;&lt;br /&gt;
Indeed, from this diagram, not only can one learn that the persistence threshold should be &lt;em&gt;~6.E20&lt;/em&gt; (from the Y-axis), but one can easily see that the nature of the distribution is different for background values below or above &quot;X=~4E20&quot; on the X-axis. This feature reflects the fact that in some regions of this data-set, the signal was below the instrument detection threshold and the result is complex noise generated by the detectors. Meaningful features may stand out from these regions, so it does not affect the choice of the persistence threshold. However, this tells us that any region with a value &lt;em&gt;&amp;lt;4E20&lt;/em&gt; cannot be considered signal and identified filaments should therefore be trimmed wherever the value is below this threshold (this can be achieved by running mse with a persitencce threhold of &lt;em&gt;6.E20&lt;/em&gt; and then running the following command on the output skeleton file: &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#trim&quot;&gt;skelconv filaments.NDskl -trimBelow 4.E20&lt;/a&gt;&lt;/em&gt;).
&lt;br /&gt;&lt;br /&gt;
Let us now come back to  our example, where no such trimming is necessary. Selecting a persistence threshold of &lt;em&gt;~4.2E21&lt;/em&gt; as discussed previously, we obtain the following filamentary structures :
&lt;br /&gt;&lt;br /&gt;
&lt;img src=&quot;http://localhost/dotclear/public/M16_skel.png&quot; alt=&quot;M16_skel.png&quot; style=&quot;display:block; margin:0 auto;&quot; /&gt;
&lt;br /&gt;&lt;br /&gt;
Remember that DisPerSE identifies structures as features of the &lt;a href=&quot;http://localhost/dotclear/index.php?post/definitions&quot;&gt;Morse-Smale complex&lt;/a&gt;. The output of &lt;em&gt;mse&lt;/em&gt; is &lt;em&gt;always&lt;/em&gt; a valid (subset of the) Morse-Smale complex. In particular, filaments are identified as the set of &lt;a href=&quot;http://localhost/dotclear/index.php?post/definitions&quot;&gt;arcs&lt;/a&gt; joining maxima and critical points, which means that pieces of filament always link a maximum and a critical point: they cannot stop in any other point, and whenever two filaments merge in points that are not maxima (i.e. bifurcation point), they remain distinct filaments until they reach a maximum (i.e. although they are not &lt;em&gt;visually&lt;/em&gt; distinct as they are infinitely close, they are still stored as separate arcs in the output file). This is illustrated on the following figure,  which shows a zoom on the identified filamentary structures with minima, saddle points and maxima represented as blue squares, green crosses and red triangles respectively:
&lt;br /&gt;&lt;br /&gt;
&lt;img src=&quot;http://localhost/dotclear/public/M16_skel_nodes.png&quot; alt=&quot;M16_skel_nodes.png&quot; style=&quot;display:block; margin:0 auto;&quot; /&gt;
&lt;br /&gt;&lt;br /&gt;
On this figure, individual pieces of filaments (i.e. arcs) link red triangles to green crosses, and several pieces of arcs may therefore overlap above bifurcation points. This is of course desirable for mathematical applications that require the computation of valid Morse-Smale complexes but can be problematic when one only needs to identify structures and measure their properties. Indeed, regions where several arcs overlap may accidentally be given a stronger weight in statistical analysis, biasing the results. One can fix this with the option &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#breakdown&quot;&gt;-breakdown&lt;/a&gt;&lt;/em&gt; of &lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv&quot;&gt;skelconv&lt;/a&gt; which is used to identify bifurcation points (adding them as dummy critical points of critical index &lt;em&gt;D+1&lt;/em&gt;) and break filaments down into individual pieces of arcs.
&lt;br /&gt;&lt;br /&gt;
Another problem is the resolution of the filaments, which is roughly equal to that of the underlying map: filaments are described as sets of small segments linking neighbor pixels together in the oupout files of &lt;em&gt;mse&lt;/em&gt;. As a result, they appear &lt;em&gt;jagged&lt;/em&gt; on the previous image and they can locally only take a discrete set of orientation. Using option &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#smooth&quot;&gt;-smooth N&lt;/a&gt;&lt;/em&gt; of &lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv&quot;&gt;skelconv&lt;/a&gt;, one can fix this by ensuring their &lt;em&gt;smoothness&lt;/em&gt; over a size of &lt;em&gt;~N pixels&lt;/em&gt;. Note that only individual arcs or pieces of arcs are smoothed when using this option, critical points and bifurcation points remain fixed. As a result, the order of options &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#smooth&quot;&gt;-smooth&lt;/a&gt;&lt;/em&gt; and &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#breakdown&quot;&gt;-breakdown&lt;/a&gt;&lt;/em&gt; on the command line is important.
&lt;br /&gt;
We choose to first smooth over  ~6 pixels and then break the network down :
&lt;br /&gt;&lt;br /&gt;&lt;/p&gt;
&lt;div style=&quot;text-align: center; border: 1px dotted gray;&quot;&gt;skelconv m16_zo_no70um_density.fits_c4.2e+21.up.NDskl -smooth 6 -breakdown&lt;/div&gt;


&lt;p&gt;&lt;br /&gt;&lt;br /&gt;
The result is illustrated on the following figure:
&lt;br /&gt;&lt;br /&gt;
&lt;img src=&quot;http://localhost/dotclear/public/M16_skelBS_nodes.png&quot; alt=&quot;M16_skelBS_nodes.png&quot; style=&quot;display:block; margin:0 auto;&quot; /&gt;
&lt;br /&gt;&lt;br /&gt;
where filaments are now smooth and bifurcation points are now identified as yellow discs : pieces of filament now link critical or bifurcation points and therefore do not overlap anymore.
&lt;br /&gt;&lt;br /&gt;
Although the resulting filaments are perfectly usable as such for analysis, other useful post-treatment are available in DisPerSE. We give in the following an example of their possible usage, but note that it should of course be adapted to the particular results one wants to obtain ...
&lt;br /&gt;&lt;br /&gt;
Options &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#trim&quot;&gt;-trimAbove&lt;/a&gt;&lt;/em&gt; and &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#trim&quot;&gt;-trimBelow&lt;/a&gt;&lt;/em&gt;  for instance allow trimming portion of arcs so that they do not have to start/stop at critical or bifurcation points (dummy critical nodes are added at their extremities). By default, these options will trim the filaments above or below a given value of the field (see e.g. the example in the &lt;a href=&quot;http://localhost/dotclear/index.php?post/regular-grid-%3A-filaments#trim_exp&quot;&gt;section&lt;/a&gt; on the persistence diagram above), but filaments can also be trimmed in any value tagged within the &lt;a href=&quot;http://localhost/dotclear/index.php?post/skeleton-formats&quot;&gt;skeleton file&lt;/a&gt; (the list of the tags can be displayed by running &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#info&quot;&gt;skelconv filaments.NDskl -info&lt;/a&gt;&lt;/em&gt; and arbitrary tags can be added with option &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#addfield&quot;&gt;-addField&lt;/a&gt;&lt;/em&gt;).
&lt;br /&gt;&lt;br /&gt;&lt;a name=&quot;robustness&quot;&gt;&lt;/a&gt;
A useful trimming function is &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/definitions#robustness&quot;&gt;robustness&lt;/a&gt;&lt;/em&gt; or &lt;em&gt;robustness_ratio&lt;/em&gt;, which is not tagged by default in the &lt;a href=&quot;http://localhost/dotclear/index.php?post/skeleton-formats&quot;&gt;skeleton files&lt;/a&gt; but can computed using option &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#robustness&quot;&gt;-robustness&lt;/a&gt;&lt;/em&gt; of &lt;em&gt;mse&lt;/em&gt; (you can use option &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#loadmsc&quot;&gt;-loadMSC&lt;/a&gt;&lt;/em&gt; to avoid recomputing the whole Morse-Smale complex):
&lt;br /&gt;&lt;br /&gt;&lt;/p&gt;
&lt;div style=&quot;text-align: center; border: 1px dotted gray;&quot;&gt;mse m16_zo_no70um_density.fits -mask M16_mask.png -cut 4.2E21 -upSkl -loadMSC m16_zo_no70um_density.fits.MSC&lt;/div&gt;


&lt;p&gt;&lt;br /&gt;&lt;br /&gt;
Robustness implements an improved version of the concept of &lt;em&gt;separatrix persistence&lt;/em&gt; (Weinkauf, T. and Gunther, D., 2009) that allows the extension of regular &lt;em&gt;persistence&lt;/em&gt; to each points of the filaments (as opposed to critical points pairs). Running the following command:
&lt;br /&gt;&lt;br /&gt;&lt;/p&gt;
&lt;div style=&quot;text-align: center; border: 1px dotted gray;&quot;&gt;skelconv m16_zo_no70um_density.fits_c4.2e+21.up.NDskl -breakdown -smooth 6 -trimBelow robustness 8E21&lt;/div&gt;


&lt;p&gt;&lt;br /&gt;&lt;br /&gt;
one can select the most meaningful subset of the filaments (you can try several robustness thresholds to see the result, but the selected threshold should probably be greater than the persistence threshold). Note that the concepts of persistence and robustness are different, and robustness basically allows selecting those subsets of persistent (topological) filaments that are most prominent. For that reason, it may make sense when trimming on robustness criteria to lower the persistence threshold (in &lt;em&gt;mse&lt;/em&gt;) so that the significant pieces of less persistent filaments can still be conserved (you should probably experiment in order to find the best compromise between persistence and robustness thresholds for a given type of data set).
In the present case, choosing a robustness threshold of 8E21 (i.e. slightly higher than the persistence threshold), one obtains the following filaments:
&lt;br /&gt;&lt;br /&gt;
&lt;img src=&quot;http://localhost/dotclear/public/M16_skelBST_inc.png&quot; alt=&quot;M16_skelBST_inc.png&quot; style=&quot;display:block; margin:0 auto;&quot; /&gt;
&lt;br /&gt;&lt;br /&gt;&lt;/p&gt;


&lt;p&gt;Another interesting post-treatment is the option &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#assemble&quot;&gt;-assemble&lt;/a&gt;&lt;/em&gt; that allows merging individual aligned pieces of filaments into larger structures. On the following picture, each piece of filament is represented with a distinct color:
&lt;br /&gt;&lt;br /&gt;
&lt;img src=&quot;http://localhost/dotclear/public/M16_skelBST_nodes_fil.png&quot; alt=&quot;M16_skelBST_nodes_fil.png&quot; style=&quot;display:block; margin:0 auto;&quot; /&gt;
&lt;br /&gt;&lt;br /&gt;
Running the following command (note once again that the order of arguments is important, try changing it ...):
&lt;br /&gt;&lt;br /&gt;&lt;/p&gt;
&lt;div style=&quot;text-align: center; border: 1px dotted gray;&quot;&gt;skelconv m16_zo_no70um_density.fits_c4.2e+21.up.NDskl -breakdown -smooth 6 -trimBelow robustness 1.0E22 -assemble 70&lt;/div&gt;


&lt;p&gt;&lt;br /&gt;&lt;br /&gt;
one can assemble pieces of filaments that form an angle smaller than &lt;em&gt;70&lt;/em&gt; degrees, resulting in the following structures, where individual pieces are straight but as long as possible:
&lt;br /&gt;&lt;br /&gt;
&lt;img src=&quot;http://localhost/dotclear/public/M16_skelBSTA_nodes_fil.png&quot; alt=&quot;M16_skelBSTA_nodes_fil.png&quot; style=&quot;display:block; margin:0 auto;&quot; /&gt;
&lt;br /&gt;&lt;br /&gt;&lt;/p&gt;</content>
    
    

    
      </entry>
  
</feed>