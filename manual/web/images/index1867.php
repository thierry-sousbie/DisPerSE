<?xml version="1.0" encoding="utf-8"?><feed xmlns="http://www.w3.org/2005/Atom"
  xmlns:dc="http://purl.org/dc/elements/1.1/"
  xmlns:wfw="http://wellformedweb.org/CommentAPI/"
  xml:lang="en">
  
  <title type="html">DisPerSE - persistent structures identification - Terminology</title>
  <subtitle type="html">Automatic identification of persistent structures in 2D or 3D.
keywords: Morse complex, topology, peak, void, source, wall, filament, cosmic web, cosmology, structure identification.</subtitle>
  <link href="http://localhost/dotclear/index.php?feed/category/Terminology/atom" rel="self" type="application/atom+xml"/>
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
    <title>Definitions</title>
    <link href="http://localhost/dotclear/index.php?post/definitions" rel="alternate" type="text/html"
    title="Definitions" />
    <id>urn:md5:8282f6bbd1d3697a7fda5befc5327eda</id>
    <published>2170-05-10T07:42:00+01:00</published>
    <updated>2012-07-16T15:25:40+02:00</updated>
    <author><name>thierry sousbie</name></author>
        <dc:subject>Terminology</dc:subject>
            
    <content type="html">    &lt;p&gt;In DisPerSE, the different types of cells that compose the &lt;a href=&quot;http://localhost/dotclear/index.php?post/general-concepts&quot;&gt;Morse-Smale complex&lt;/a&gt; are used to identify different types of structures. In this document, we adopt the following terminology:
&lt;br /&gt;
&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;- &lt;ins&gt;A &lt;strong&gt;critical point&lt;/strong&gt;&lt;/ins&gt;&lt;br /&gt; is a point where the gradient is null. In N dimensional space, there are N types of critical point, characterized by a critical index.&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;- &lt;ins&gt;A &lt;strong&gt;maximum&lt;/strong&gt;&lt;/ins&gt;&lt;br /&gt; is a critical point of critical index N in an N dimensional space.&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;- &lt;ins&gt;A &lt;strong&gt;minimum&lt;/strong&gt;&lt;/ins&gt;&lt;br /&gt; is a critical point of critical index 0.&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;- &lt;ins&gt;A &lt;strong&gt;k-saddle point&lt;/strong&gt;&lt;/ins&gt;&lt;br /&gt; is a critical point of critical index k that is neither a minimum nor a maximum.&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;- &lt;ins&gt;A &lt;strong&gt;d-simplex&lt;/strong&gt;&lt;/ins&gt;&lt;br /&gt; is the most simple geometrical figure in &lt;em&gt;d&lt;/em&gt; dimensions and consists of &lt;em&gt;d+1&lt;/em&gt; points linked together (i.e a point, segment, triangle, tetrahedron ... for d=0,1,2,3,...). In DisPerSE, space is decomposed into d-simplices and each critical point is associated to a given &lt;em&gt;d&lt;/em&gt;-simplex (by default, maxima are 0-simplices, saddle-points are 1-simplices, ...).&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;- &lt;ins&gt;A &lt;strong&gt;persistence pair&lt;/strong&gt;&lt;/ins&gt;&lt;br /&gt; is a pair of critical points with critical index difference of 1. A persistence pair represents a topological component of the function. More details &lt;a href=&quot;http://localhost/dotclear/index.php?post/Persistence-and-simplification&quot;&gt;here&lt;/a&gt;.&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;- &lt;ins&gt; &lt;strong&gt;Persistence&lt;/strong&gt;&lt;/ins&gt;&lt;br /&gt; is defined as the difference (or ratio) of the value at the two critical points in a persistence pair. Persistence is always positive and represents the importance of the topological feature defined by the persistence pair : low persistence pairs are sensible to changes in the value of the function, and the corresponding critical points are easily destroyed, even by small perturbations, while high persistence pairs are robust. Removing low persistence pairs from the &lt;a href=&quot;http://localhost/dotclear/index.php?post/general-concepts&quot;&gt;Morse-Smale complex&lt;/a&gt; is a good way to filter noise and/or non-meaningfull structures. More details &lt;a href=&quot;http://localhost/dotclear/index.php?post/Persistence-and-simplification&quot;&gt;here&lt;/a&gt;.&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;robustness&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;- &lt;ins&gt;&lt;strong&gt;Robustness&lt;/strong&gt;&lt;/ins&gt;&lt;br /&gt; is a local measure of how contrasted the critical points and filaments are with respect to their background. In a sense, robustness is a continuous version of persistence which extends to filaments (persistence applies to pairs of critical points only). What we call robustness here is actually an improved version of the &lt;em&gt;separatrix persistence&lt;/em&gt; described in &lt;em&gt;Weinkauf, T. and Gunther, D., 2009&lt;/em&gt;. See options &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#trim&quot;&gt;-trimBelow&lt;/a&gt;&lt;/em&gt; of &lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv&quot;&gt;skelconv&lt;/a&gt;, &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#robustness&quot;&gt;-robustness&lt;/a&gt;&lt;/em&gt; of &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse&quot;&gt;mse&lt;/a&gt;&lt;/em&gt; and the &lt;a href=&quot;http://localhost/dotclear/index.php?post/regular-grid-%3A-filaments#robustness&quot;&gt;tutorial section&lt;/a&gt; for more details and possible applications.&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;- &lt;ins&gt;An &lt;strong&gt;integral line&lt;/strong&gt;&lt;/ins&gt;&lt;br /&gt; is a curve tangent to the gradient field in every point. There exist exactly one integral line going through every non critical point of the domain of definition, and gradient lines must start and end at critical points (i.e. where the gradient is null).&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;- &lt;ins&gt;A &lt;strong&gt;descending k-manifold&lt;/strong&gt;&lt;/ins&gt;&lt;br /&gt; is a k dimensional region of space defined as the set of points from which following minus the gradient leads to the same critical point of critical index k, with N the dimension of the space (e.g. a descending N-manifold originates from a maximum). More details &lt;a href=&quot;http://localhost/dotclear/index.php?post/general-concepts&quot;&gt;here&lt;/a&gt;.&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;- &lt;ins&gt;An &lt;strong&gt;ascending k-manifold&lt;/strong&gt;&lt;/ins&gt;&lt;br /&gt; is a k dimensional region of space defined as the set of points from which following the gradient leads to the same critical point of critical index N-k, where N is the dimension of the space (e.g. an ascending N-manifold originates from a minimum). More details &lt;a href=&quot;http://localhost/dotclear/index.php?post/general-concepts&quot;&gt;here&lt;/a&gt;.&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;- &lt;ins&gt;A &lt;strong&gt;Morse p-cell&lt;/strong&gt;&lt;/ins&gt;&lt;br /&gt; is a p-dimensional region of space defined by the intersection of an ascending and a descending manifold. All Integral lines in a morse p-cell originate and lead to the same critical points. More details &lt;a href=&quot;http://localhost/dotclear/index.php?post/general-concepts&quot;&gt;here&lt;/a&gt;.&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;- &lt;ins&gt;A &lt;strong&gt;Morse-Smale complex&lt;/strong&gt;&lt;/ins&gt;&lt;br /&gt; of a Morse function is the tessellation of space induced by its gradient : it divides space into regions where the integral lines originate from and lead to the same critical points. Its elements are the Morse p-cells. More details &lt;a href=&quot;http://localhost/dotclear/index.php?post/general-concepts&quot;&gt;here&lt;/a&gt;.&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;- &lt;ins&gt;An &lt;strong&gt;arc&lt;/strong&gt;&lt;/ins&gt;&lt;br /&gt; is a Morse 1-cell : a 1D curve traced by the integral line that joins two critical points together. There is always exactly two arcs originating from a saddle point and leading to an extremum (i.e. a minimum or a maximum). In 3D and more however, the number of arcs linking two saddle point together is arbitrary. More details &lt;a href=&quot;http://localhost/dotclear/index.php?post/general-concepts&quot;&gt;here&lt;/a&gt;.&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;- &lt;ins&gt;A &lt;strong&gt;skeleton&lt;/strong&gt;&lt;/ins&gt;&lt;br /&gt; is a set of arcs. We call &lt;em&gt;upper-&lt;/em&gt;, &lt;em&gt;lower-&lt;/em&gt;, and &lt;em&gt;inter-&lt;/em&gt; skeleton the set of arcs that lead or originate from exactly one maximum, exactly one minimum, or none of the two respectively.&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;- &lt;ins&gt;A &lt;strong&gt;filament&lt;/strong&gt;&lt;/ins&gt;&lt;br /&gt; is a 1-dimensional structure : an ascending or descending 1-manifold. A filament consists of the set of two arcs originating from a given saddle point and joining two extrema together (i.e. the yellow+green arcs on lower right frame of &lt;a href=&quot;http://localhost/dotclear/index.php?post/general-concepts&quot;&gt;figure 1&lt;/a&gt;). Unless explicitly stated, filaments join maxima together and we sometime call void filaments or anti-filaments the filaments that link minima together.&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;- &lt;ins&gt;A &lt;strong&gt;wall&lt;/strong&gt;&lt;/ins&gt;&lt;br /&gt; is a 2-dimensional structure : an ascending or descending 2-manifold. Walls exist only in 3D and more and delimit the voids or peak patches.&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;- &lt;ins&gt;A &lt;strong&gt;void&lt;/strong&gt;&lt;/ins&gt;&lt;br /&gt; is a N dimensional structure (&lt;em&gt;N&lt;/em&gt; being the dimension of space): an ascending N-manifold that originates from a given minimum of the field. A void is the set of points from which following the gradient leads to a given minimum.&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;- &lt;ins&gt;A &lt;strong&gt;peak patch&lt;/strong&gt;&lt;/ins&gt;&lt;br /&gt; is a N dimensional structure (&lt;em&gt;N&lt;/em&gt; being the dimension of space): a descending N-manifold that originates from a given maximum of the field. Peak patches are the dual of voids: the set of points from which following the opposite of the gradient leads to the same maxima.&lt;/li&gt;
&lt;/ul&gt;</content>
    
    

    
      </entry>
  
</feed>