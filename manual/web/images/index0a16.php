<?xml version="1.0" encoding="utf-8"?><feed xmlns="http://www.w3.org/2005/Atom"
  xmlns:dc="http://purl.org/dc/elements/1.1/"
  xmlns:wfw="http://wellformedweb.org/CommentAPI/"
  xml:lang="en">
  
  <title type="html">DisPerSE - persistent structures identification - Overview</title>
  <subtitle type="html">Automatic identification of persistent structures in 2D or 3D.
keywords: Morse complex, topology, peak, void, source, wall, filament, cosmic web, cosmology, structure identification.</subtitle>
  <link href="http://localhost/dotclear/index.php?feed/category/Overview/atom" rel="self" type="application/atom+xml"/>
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
    <title>Persistence and simplification</title>
    <link href="http://localhost/dotclear/index.php?post/Persistence-and-simplification" rel="alternate" type="text/html"
    title="Persistence and simplification" />
    <id>urn:md5:66c05bde8c3be493a1b4c025fc072c47</id>
    <published>2160-06-10T07:41:00+01:00</published>
    <updated>2012-06-07T14:36:10+02:00</updated>
    <author><name>thierry sousbie</name></author>
        <dc:subject>Overview</dc:subject>
            
    <content type="html">    &lt;p&gt;Persistence itself is a relatively simple but powerful concept. To study the topology of a function, one can measure how the topology of its excursion sets (i.e. the set of points with value higher than a given threshold) evolves when the threshold is continuously and monotonically changing. Whenever the threshold crosses the value of a critical point, the topology of the excursion change. Supposing that the threshold is sweeping the values of a 1D function from high to low, whenever it crosses the value of a maximum, a new component appears in the excursion, while two components merge (i.e. one is destroyed) whenever the threshold crosses the value of a minimum. This concept can be extended to higher dimensions (i.e. creation/destruction of hole, spherical shells, ....) and in general, whenever a topological component is created at a critical point, the critical point is labeled positive, while it is labeled negative if it destroys a topological component. Using this definition, topological components of a function can be represented by pairs of positive and negative critical points called persistence pairs. The absolute difference of the value of the critical points in a pair is called its persistence : it represents the lifetime of the corresponding topological component within the excursion set.
&lt;br /&gt;
&lt;br /&gt;&lt;/p&gt;
&lt;figure style=&quot;text-align:center;font-family:times;color:rgb(32,32,32);font-size:0.9em;&quot;&gt;
Figure 1: Persistence pairs&lt;/br&gt;
&lt;img src=&quot;http://localhost/dotclear/public/persistence1D.png&quot; alt=&quot;Blah&quot;/&gt;&lt;br&gt;
&lt;legend style=&quot;text-align:justify;font-family:times;color:rgb(32,32,32);font-size:1em;line-height:90%;margin:0px 20px 0px 20px;&quot;&gt;
&lt;i style=&quot;font-size:0.9em;&quot;&gt;Upper part:&lt;/i&gt; Two functions with identical topology (i.e. two peaks).
&lt;i style=&quot;font-size:0.9em;&quot;&gt;Lower part:&lt;/i&gt; Changes in topology of excursion sets (sets of point with value lower than a decreasing threshold, see main text above). In 1D, new components are created at maxima and destroyed at minima. Two persistence pairs represented by green lines account for the two peaks in each function. The length of the green lines corresponds to the persistence of the pairs, which correctly accounts for the the fact that A has a small bump on the top of a peak (one high and one low persistence pair), while B has two peaks (two high persistence pairs).
&lt;/legend&gt;
&lt;/figure&gt;


&lt;p&gt;&lt;br /&gt;
&lt;br /&gt;
The concept of persistence is powerful because it yields a simple way to measure how robust topological components are to local modifications of a function values. Indeed, noise can only affect a function's topology by creating or destroying topological components of persistence lower that its local amplitude. Therefore, it suffice to know the amplitude of noise to decide which components certainly belong to an underlying function and which may have been affected (i.e. created or destroyed) by noise. In DisPerSE, a persistence threshold can be specified (see options &lt;em&gt;-nsig&lt;/em&gt; and &lt;em&gt;-cut&lt;/em&gt; of &lt;a href=&quot;http://localhost/dotclear/index.php?post/mse&quot;&gt;mse&lt;/a&gt;) to remove topological components with persistence lower than the threshold and therefore filter noise from the Morse-Smale complex (see figure 2 below).
&lt;br /&gt;
&lt;br /&gt;&lt;/p&gt;
&lt;figure style=&quot;text-align:center;font-family:times;color:rgb(32,32,32);font-size:0.9em;&quot;&gt;
Figure 1: Morse-smale complexsimplification&lt;/br&gt;
&lt;img src=&quot;http://localhost/dotclear/public/persistence_simp2D.jpg&quot; alt=&quot;Blah&quot;/&gt;&lt;br&gt;
&lt;legend style=&quot;text-align:justify;font-family:times;color:rgb(32,32,32);font-size:1em;line-height:90%;margin:0px 20px 0px 20px;&quot;&gt;
Simplification of a low persistence pair (formed by maximum 19 and saddle point 18). Disks represent critical points (minima, saddle-points and maxima in blue, green and red respectively) and black lines show the arcs of the Morse-Smale complex.
&lt;/legend&gt;
&lt;/figure&gt;


&lt;p&gt;&lt;br /&gt;
&lt;br /&gt;
A very useful way to set the persistence threshold is to plot a persistence diagram, in which all the persistence pairs are represented by points with coordinates the value at the critical points in the pair (see option &lt;em&gt;-interactive&lt;/em&gt; of &lt;a href=&quot;http://localhost/dotclear/index.php?post/mse&quot;&gt;mse&lt;/a&gt; or see &lt;a href=&quot;http://localhost/dotclear/index.php?post/pdview&quot;&gt;pdview&lt;/a&gt; and read the &lt;em&gt;tutorial&lt;/em&gt; section to learn how to compute and use persistence diagrams).&lt;/p&gt;</content>
    
    

    
      </entry>
    
  <entry>
    <title>General concepts of Morse theory</title>
    <link href="http://localhost/dotclear/index.php?post/general-concepts" rel="alternate" type="text/html"
    title="General concepts of Morse theory" />
    <id>urn:md5:ed25efa8ed598b54d952044d668b4fbe</id>
    <published>2160-05-10T07:41:00+01:00</published>
    <updated>2012-06-07T14:54:37+02:00</updated>
    <author><name>thierry sousbie</name></author>
        <dc:subject>Overview</dc:subject>
            
    <content type="html">    &lt;p&gt;In DisPerSE, structures are identified as components of the Morse-Smale complex of an input function defined over a - possibly bounded - manifold. The Morse-Smale complex of a real valued so-called Morse function is a construction of Morse theory which captures the relationship between the gradient of the function, its topology, and the topology of the manifold it is defined over.
&lt;br /&gt;
&lt;br /&gt;&lt;/p&gt;
&lt;figure style=&quot;text-align:center;font-family:times;color:rgb(32,32,32);font-size:0.9em;&quot;&gt;
Figure 1: the Morse-Smale complex&lt;/br&gt;
&lt;img src=&quot;http://localhost/dotclear/public/morse_illustration.jpg&quot; alt=&quot;Blah&quot;/&gt;&lt;br&gt;
&lt;legend style=&quot;text-align:justify;font-family:times;color:rgb(32,32,32);font-size:1em;line-height:90%;margin:0px 20px 0px 20px;&quot;&gt;
&lt;i style=&quot;font-size:0.9em;&quot;&gt;Upper left:&lt;/i&gt; critical points (minima, saddle points and maxima pictured as blue, green and red disks), and three integral lines (pink curves) of a Morse function. Black arrows show the gradient of that function.
&lt;i style=&quot;font-size:0.9em;&quot;&gt;Upper right:&lt;/i&gt; ascending 2-manifolds : the set of points belonging to integral lines whose destination is the same minimum (critical point of index 0).
&lt;i style=&quot;font-size:0.9em;&quot;&gt;Lower left:&lt;/i&gt; descending 2-manifolds : the set of points belonging to integral lines whose origin is the same maximum (critical point of index 2).
&lt;i style=&quot;font-size:0.9em;&quot;&gt;Lower right:&lt;/i&gt; the Morse-Smale complex : a natural tesselation of space into cells induced by the gradient fo the function. Each cell is the set of points belonging to integral lines whose origin and destination are identical (i.e. each cell is the intersection of an ascending and a descending manifold). The purple region is a 2-cell: intersection of an ascending and a descending 2-manifold (red and blue regions) where all field lines have the same orgin and destination (a minimum and a maxium). The yellow curve is a 1-cell (also called an arc): the intersection of and ascending 2-manifold (blue region) and a descending 1-manifolds (green+yellow curves, originating from the same saddle point).
&lt;/legend&gt;
&lt;/figure&gt;


&lt;p&gt;&lt;br /&gt;
&lt;br /&gt;
Two central notions in Morse theory are that of &lt;em&gt;critical point&lt;/em&gt; and &lt;em&gt;integral line&lt;/em&gt; (also called &lt;em&gt;field line&lt;/em&gt;) : (see the upper left frame of figure 1)
&lt;br /&gt;
&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;Critical points&lt;/strong&gt; are the discrete set of points where the gradient of the function is null. For a function defined over a 2D space, there are three types of critical points (4 in 3D, ...), classified by their critical index. In 2D, minima have a critical index of 0, saddle points have a critical index of 1 and maxima have a critical index of 2. In 3D and more, different types of saddle points exist, one for each non extremal critical index.&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;Integral lines&lt;/strong&gt; are curves tangent to the gradient field in every point. There exist exactly one integral line going through every non critical point of the domain of definition, and gradient lines must start and end at critical points (i.e. where the gradient is null).&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;
Because integral lines cover all space (there is exactly one critical line going through every point of space) and their extremities are critical points, they induce a tessellation of space into regions called ascending (resp. descending) &lt;em&gt;k-manifolds&lt;/em&gt; where all the field lines originate (respectively lead) from the same critical point (see ascending and descending 2-manifolds on figure 1, upper right and lower left panels). The number of dimensions &lt;em&gt;k&lt;/em&gt; of the regions spanned by a k-manifold depend directly on the critical index of the corresponding critical point: descending k-manifolds originate from critical points of critical index &lt;em&gt;k&lt;/em&gt; while the critical index is &lt;em&gt;N-k&lt;/em&gt; for ascending manifolds, with &lt;em&gt;N&lt;/em&gt; the dimension of space.
&lt;br /&gt;
&lt;br /&gt;
The set of all ascending (or descending) manifolds is called the &lt;em&gt;Morse complex&lt;/em&gt; of the function. The &lt;em&gt;Morse-Smale complex&lt;/em&gt; is an extension of this concept: the tessellation of space into regions called &lt;em&gt;p-cells&lt;/em&gt; where all the integral lines have the same origin and destination (see figure 1, lower right frame). Each p-cell of the Morse-Smale complex is the intersection of an ascending and a descending manifold and the Morse-Smale complex itself is a natural tessellation of space induced by the gradient of the function. Figure 2 below illustrates how components of the Morse-Smale complex can be used to identify structures in a 3D distribution.
&lt;br /&gt;
&lt;br /&gt;&lt;/p&gt;
&lt;figure style=&quot;text-align:center;font-family:times;color:rgb(32,32,32);font-size:0.9em;&quot;&gt;
Figure 2: 3D structures identified as component of the Morse-Smale complex
&lt;img src=&quot;http://localhost/dotclear/public/structures_as_manifolds.jpg&quot; alt=&quot;Blah&quot;/&gt;&lt;br&gt;
&lt;legend style=&quot;text-align:justify;font-family:times;color:rgb(32,32,32);font-size:1em;line-height:90%;margin:0px 20px 0px 20px;&quot;&gt;
&lt;i style=&quot;font-size:0.9em;&quot;&gt;Upper left:&lt;/i&gt; Density distribution of dark matter in a chunk of the Universe represented by tracer particles from a N-Body cosmological simulation.
&lt;i style=&quot;font-size:0.9em;&quot;&gt;Upper right:&lt;/i&gt; Ascending 3-manifolds tracing the voids
&lt;i style=&quot;font-size:0.9em;&quot;&gt;Lower left:&lt;/i&gt; Ascending 2-manifolds tracing the walls
&lt;i style=&quot;font-size:0.9em;&quot;&gt;Lower right:&lt;/i&gt; The set of arcs with one maximum at their extremity, also called upper skeleton, tracing the filamentary structures. The maxima, not represented here, identify dark matter halos onto which filaments plug.
&lt;/legend&gt;
&lt;/figure&gt;


&lt;p&gt;&lt;br /&gt;
&lt;br /&gt;&lt;/p&gt;</content>
    
    

    
      </entry>
    
  <entry>
    <title>Examples of applications</title>
    <link href="http://localhost/dotclear/index.php?post/Examples-of-applications" rel="alternate" type="text/html"
    title="Examples of applications" />
    <id>urn:md5:d4284311032bb57e5c1192fa135a66bb</id>
    <published>2155-05-11T09:32:00+01:00</published>
    <updated>2013-01-16T02:47:16+01:00</updated>
    <author><name>thierry sousbie</name></author>
        <dc:subject>Overview</dc:subject>
            
    <content type="html">    &lt;p&gt;Here are a few examples of what DisPerSE can do (hover the cursor over the pictures for a short description) :
&lt;br /&gt;
&lt;br /&gt;
&lt;img src=&quot;http://localhost/dotclear/public/.simu2D_filter2-B_m.jpg&quot; alt=&quot;cosmic web&quot; style=&quot;float:left; margin: 0 1em 1em 0;&quot; title=&quot;Slice of a discretely sampled dark matter distribution in a numerical simulation of a 50Mpc chunk of the Universe. The filamentary structure of the cosmic web is clearly visible.&quot; /&gt;
&lt;img src=&quot;http://localhost/dotclear/public/.simu2D_filter2-A5_m.jpg&quot; alt=&quot;filamentary structures&quot; style=&quot;float:right; margin: 0 0 1em 1em;&quot; title=&quot;Identified filamentary structure and corresponding critical points (maxima in red, minima in blue and saddle points in yellow) &quot; /&gt;
&lt;br /&gt;
&lt;br /&gt;
&lt;img src=&quot;http://localhost/dotclear/public/.simu2D_filter-B_m.jpg&quot; alt=&quot;zoom&quot; style=&quot;float:left; margin: 0 1em 1em 0;&quot; title=&quot;Zoom on a clump, the color corresponds to the density computed by DTFE&quot; /&gt;
&lt;img src=&quot;http://localhost/dotclear/public/.simu2D_filter-A5_m.jpg&quot; alt=&quot;filaments connecting on a dark matter halo&quot; style=&quot;float:right; margin: 0 0 1em 1em;&quot; title=&quot;Persistent filaments connecting to a clump, note how maxima (red triangles) correctly identify bound structures&quot; /&gt;
&lt;br /&gt;
&lt;br /&gt;
&lt;img src=&quot;http://localhost/dotclear/public/CMB_skl_small.png&quot; alt=&quot;CMB filaments&quot; style=&quot;display:block; margin:0 auto;&quot; title=&quot; Filamentary structures in the cosmic microwave background (CMB). Structures were directly identified over the Healpix tessellation of the sphere&quot; /&gt;
&lt;br /&gt;
&lt;br /&gt;
&lt;img src=&quot;http://localhost/dotclear/public/.F3_simu250_large_simu_void_skel_crit_m.jpg&quot; alt=&quot;&quot; style=&quot;float:left; margin: 0 1em 1em 0;&quot; title=&quot;Filamentary structures and a void (lower right) identified in the simulated distribution of dark matter. Structures where identified over the full 3D distribution, but only a 20Mpc thick slice is shown for convenience.&quot; /&gt;
&lt;img src=&quot;http://localhost/dotclear/public/.F3_simu250Z_void_simu_m.jpg&quot; alt=&quot;&quot; style=&quot;float:left; margin: 0 1em 1em 0;&quot; title=&quot;Zoom on the void (i.e. low density region). Color corresponds to the density on the surface of the void.&quot; /&gt;
&lt;br /&gt;
&lt;br /&gt;
&lt;img src=&quot;http://localhost/dotclear/public/herschel_fil.png&quot; alt=&quot;herschel_fil.png&quot; style=&quot;display:block; margin:0 auto;&quot; title=&quot; Star forming regions in the interstellar medium as observed by HERSCHEL.&quot; /&gt;&lt;/p&gt;</content>
    
    

    
      </entry>
    
  <entry>
    <title>Purpose</title>
    <link href="http://localhost/dotclear/index.php?post/General-overview" rel="alternate" type="text/html"
    title="Purpose" />
    <id>urn:md5:1175fd3311f0db0202aae88a8c8a7288</id>
    <published>2150-05-09T03:38:00+01:00</published>
    <updated>2012-10-02T16:22:54+02:00</updated>
    <author><name>thierry sousbie</name></author>
        <dc:subject>Overview</dc:subject>
            
    <content type="html">    &lt;p&gt;DisPerSE stands for &lt;q&gt;Discrete Persistent Structures Extractor&lt;/q&gt; and its main purpose is to identify persistent topological features such as peaks, voids, walls and in particular filamentary structures within sampled distributions in 2D, 3D, and possibly more ...
&lt;br /&gt;
&lt;br /&gt;
Although it was initially developed with cosmology in mind (for the study of the properties of filamentary structures in the so called comic web of galaxy distribution over large scales in the Universe), the present version is quite versatile and should be useful for any application where a robust structure identification is required, for segmentation or for studying the topology of sampled functions (like computing persistent Betti numbers for instance).
&lt;br /&gt;
&lt;br /&gt;
DisPerSE is able to deal directly with noisy datasets using the concept of persistence (a measure of the robustness of topological features) and can work indifferently on many kinds of cell complex (such as structured and unstructured grids, 2D manifolds embedded within a 3D space, discrete point samples using delaunay tesselation, Healpix tesselations of the sphere, ...). The only constraint is that the distribution must be defined over a manifold, possibly with boundaries.
&lt;br /&gt;
&lt;br /&gt;&lt;/p&gt;</content>
    
    

    
      </entry>
  
</feed>