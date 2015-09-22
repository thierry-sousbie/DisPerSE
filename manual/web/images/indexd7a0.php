<?xml version="1.0" encoding="utf-8"?><feed xmlns="http://www.w3.org/2005/Atom"
  xmlns:dc="http://purl.org/dc/elements/1.1/"
  xmlns:wfw="http://wellformedweb.org/CommentAPI/"
  xml:lang="en">
  
  <title type="html">DisPerSE - persistent structures identification - Manual</title>
  <subtitle type="html">Automatic identification of persistent structures in 2D or 3D.
keywords: Morse complex, topology, peak, void, source, wall, filament, cosmic web, cosmology, structure identification.</subtitle>
  <link href="http://localhost/dotclear/index.php?feed/category/Manual/atom" rel="self" type="application/atom+xml"/>
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
    <title>fieldconv</title>
    <link href="http://localhost/dotclear/index.php?post/fieldconv" rel="alternate" type="text/html"
    title="fieldconv" />
    <id>urn:md5:f7f21449ea1a899b666bf00d747ec18d</id>
    <published>2260-05-09T06:51:00+01:00</published>
    <updated>2012-05-22T18:20:22+02:00</updated>
    <author><name>thierry sousbie</name></author>
        <dc:subject>Manual</dc:subject>
            
    <content type="html">    &lt;p&gt;&lt;strong&gt;fieldconv&lt;/strong&gt; is used to display information about &lt;a href=&quot;http://localhost/dotclear/index.php?post/field-formats&quot;&gt;field format&lt;/a&gt; files encoding regular grid or point set coordinates, and convert them to other formats. Run &lt;em&gt;fieldconv&lt;/em&gt; without argument to see a list of available input and output formats.&lt;/p&gt;


&lt;p&gt;&lt;br /&gt;
&lt;br /&gt;&lt;/p&gt;


&lt;p&gt;&lt;ins&gt;Usage:&lt;/ins&gt;&lt;/p&gt;

&lt;pre&gt; &lt;strong&gt;fieldconv&lt;/strong&gt; &amp;lt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/fieldconv#filename&quot;&gt;filename&lt;/a&gt;&amp;gt; 
           [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/fieldconv#outName&quot;&gt;outName&lt;/a&gt; &amp;lt;output filename&amp;gt;] [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/fieldconv#outDir&quot;&gt;outDir&lt;/a&gt; &amp;lt;dir&amp;gt;]
           [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/fieldconv#info&quot;&gt;info&lt;/a&gt;] [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/fieldconv#to&quot;&gt;to&lt;/a&gt; &amp;lt;format&amp;gt;]&lt;/pre&gt;

&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;filename&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;&amp;lt;filename&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;The name of a file containing a regular grid or point set coordinates in a readable &lt;a href=&quot;http://localhost/dotclear/index.php?post/field-formats&quot;&gt;field format&lt;/a&gt;.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;outName&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-outName &amp;lt;fname&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Specifies the base name of the output file (extensions are added to this base name depending on the output type).&lt;br /&gt;&lt;ins&gt;Default value&lt;/ins&gt;: the name of the input file.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;outDir&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-outDir &amp;lt;dir&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Specifies the output directory. &lt;br /&gt;&lt;ins&gt;Default value&lt;/ins&gt;: the current working directory.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;info&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-info&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Prints information on the input file.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;to&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-to &amp;lt;format&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Outputs a file in the selected writable &lt;a href=&quot;http://localhost/dotclear/index.php?post/field-formats&quot;&gt;field format&lt;/a&gt;. A list of possible parameter values can be obtained by running &lt;em&gt;fieldconv&lt;/em&gt; without any argument.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;</content>
    
    

    
      </entry>
    
  <entry>
    <title>skelconv</title>
    <link href="http://localhost/dotclear/index.php?post/skelconv" rel="alternate" type="text/html"
    title="skelconv" />
    <id>urn:md5:8a823545f404f7af8e12a4ea63bba9ce</id>
    <published>2250-05-09T06:49:00+01:00</published>
    <updated>2013-01-11T04:38:55+01:00</updated>
    <author><name>thierry sousbie</name></author>
        <dc:subject>Manual</dc:subject>
            
    <content type="html">    &lt;p&gt;&lt;strong&gt;skelconv&lt;/strong&gt;  is used to post-treat, convert and display information about &lt;a href=&quot;http://localhost/dotclear/index.php?post/skeleton-formats&quot;&gt;skeleton format&lt;/a&gt; files. Run &lt;em&gt;skelconv&lt;/em&gt; without argument to see a list of available input and output formats. Note that skeleton files can also be converted to &lt;a href=&quot;http://localhost/dotclear/index.php?post/networks&quot;&gt;NDnet format&lt;/a&gt; unstructured networks readable by &lt;a href=&quot;http://localhost/dotclear/index.php?post/netconv&quot;&gt;netconv&lt;/a&gt;.
&lt;br /&gt;&lt;ins&gt;&lt;strong&gt;nb&lt;/strong&gt;&lt;/ins&gt;: in skelconv, the order of the argument on the command line is important and certain post-treatment do not commute.  In particular, be careful when using options &lt;em&gt;-&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#breakdown&quot;&gt;breakdown&lt;/a&gt;&lt;/em&gt;, &lt;em&gt;-&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#trim&quot;&gt;trim&lt;/a&gt;&lt;/em&gt;, &lt;em&gt;-&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#smooth&quot;&gt;smooth&lt;/a&gt;&lt;/em&gt; or &lt;em&gt;-&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#assemble&quot;&gt;assemble&lt;/a&gt;&lt;/em&gt; together as for instance &lt;em&gt;-smooth 3 -breakdown&lt;/em&gt; or &lt;em&gt;-breakdown -smooth 3&lt;/em&gt; will not produce the same output.
&lt;br /&gt;
&lt;br /&gt;&lt;/p&gt;


&lt;p&gt;&lt;ins&gt;Usage:&lt;/ins&gt;&lt;/p&gt;

&lt;pre&gt; &lt;strong&gt;skelconv&lt;/strong&gt; &amp;lt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#filename&quot;&gt;filename&lt;/a&gt;&amp;gt; 
          [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#outName&quot;&gt;outName&lt;/a&gt; &amp;lt;output filename&amp;gt;] [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#outDir&quot;&gt;outDir&lt;/a&gt; &amp;lt;dir&amp;gt;] [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#noTags&quot;&gt;noTags&lt;/a&gt;]
          [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#smooth&quot;&gt;smooth&lt;/a&gt; &amp;lt;Ntimes=0&amp;gt;] [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#breakdown&quot;&gt;breakdown&lt;/a&gt;]
          [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#assemble&quot;&gt;assemble&lt;/a&gt; [[&amp;lt;field_name&amp;gt;] &amp;lt;threshold&amp;gt;] &amp;lt;angle&amp;gt;]
          [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#trim&quot;&gt;trimAbove&lt;/a&gt; [&amp;lt;field_name&amp;gt;] &amp;lt;threshold&amp;gt;] 
          [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#trim&quot;&gt;trimBelow&lt;/a&gt; [&amp;lt;field_name&amp;gt;] &amp;lt;threshold&amp;gt;] 
          [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#rmBoundary&quot;&gt;rmBoundary&lt;/a&gt;] [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#rmOutside&quot;&gt;rmOutside&lt;/a&gt;]
          [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#toRaDecZ&quot;&gt;toRaDecZ&lt;/a&gt;] [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#toRaDecDist&quot;&gt;toRaDecDist&lt;/a&gt;]
          [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#toFITS&quot;&gt;toFITS&lt;/a&gt; [&amp;lt;Xres&amp;gt;] [&amp;lt;Yres&amp;gt;] [&amp;lt;Zres&amp;gt;]]
          [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#info&quot;&gt;info&lt;/a&gt;] [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#addField&quot;&gt;addField&lt;/a&gt; &amp;lt;filename&amp;gt; &amp;lt;field_name&amp;gt;] [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#to&quot;&gt;to&lt;/a&gt; &amp;lt;format&amp;gt;]&lt;/pre&gt;

&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;filename&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;&amp;lt;filename&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;The name of a file containing a skeleton such as skeletons (filaments) output by &lt;a href=&quot;http://localhost/dotclear/index.php?post/mse&quot;&gt;mse&lt;/a&gt; in a readable &lt;a href=&quot;http://localhost/dotclear/index.php?post/skeleton-formats&quot;&gt;skeleton format&lt;/a&gt;.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;outName&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-outName &amp;lt;fname&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Specifies the base name of the output file (extensions are added to this base name depending on the output type).&lt;br /&gt;&lt;ins&gt;Default value&lt;/ins&gt;: the name of the input file.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;outDir&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-outDir &amp;lt;dir&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Specifies the output directory. &lt;br /&gt;&lt;ins&gt;Default value&lt;/ins&gt;: the current working directory.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;noTags&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-noTags&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Prevents &lt;em&gt;skelconv&lt;/em&gt; from adding trailing extensions to the output filename. Note that the last extension correponding to the file format is still added.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;smooth&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-smooth &amp;lt;Ntimes=0&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Smooth the skeleton N times. When using this option, the nodes (i.e. critical points) are fixed and the filaments are smoothed by averaging the coordinates of each points along a filament with that of its two neighbors. Smoothing &lt;em&gt;N&lt;/em&gt; times effectively make filaments smooth over &lt;em&gt;~N&lt;/em&gt; sampling points.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;breakdown&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-breakdown&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;By default, skeletons are composed of arcs linking critical points together, so an arc will always start and stop at a critical point. Different arcs with same destination may partially overlap though, between a so-called bifurcation and a critical point. If N arcs overlap in a given place, the N segments will describe their geometry in that place : this is topologically correct but may not be desirable when computing statistics along the filaments as it would artificially increase the weight of those regions. Using &lt;em&gt;-breakdown&lt;/em&gt;, bifurcation points are replaced by fake critical points with critical index D+1 (where D is the number of dimensions), and infinitely close pieces of arcs are merged. This option should most probably be used before computing any statistical quantity along arcs ...&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;assemble&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-assemble [[&amp;lt;field_name&amp;gt;] &amp;lt;threshold&amp;gt;] &amp;lt;angle&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Assemble arcs into longer filaments, with the constraint that they do not form an angle larger than &lt;em&gt;&amp;lt;angle&amp;gt;&lt;/em&gt; (expressed in degrees) and optionally do not go below threshold &lt;em&gt;&amp;lt;threshold&amp;gt;&lt;/em&gt; in the field &lt;em&gt;&amp;lt;field_name&amp;gt;&lt;/em&gt;. See option &lt;em&gt;-&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#trim&quot;&gt;trim&lt;/a&gt;&lt;/em&gt; for explanations on the first two arguments. When using this option, the algorithm will try to find the longest possible aligned arcs and will join them, removing critical points and creating only straight filaments. Option &lt;em&gt;-&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#breakdown&quot;&gt;breakdown&lt;/a&gt;&lt;/em&gt; should almost always be used before using this (e.g. &lt;em&gt;skelconv filename.NDskl -breakdown -assemble 0 45&lt;/em&gt;).&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;trim&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-trimAbove/trimBelow [&amp;lt;field_name&amp;gt;] &amp;lt;threshold&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Trims the regions of the skeleton above or below threshold &lt;em&gt;&amp;lt;threshold&amp;gt;&lt;/em&gt; and add nodes (fake critical points of index D+1) at the extremity of trimmed arcs. By default, the skeleton is trimmed according to &lt;em&gt;field_value&lt;/em&gt; (i.e. the value of the function from which the skeleton was computed). One can specify a different function with &lt;em&gt;&amp;lt;field_name&amp;gt;&lt;/em&gt; (a list of possible fields may be obtained running &lt;em&gt;skelconv filename.NDskl -info&lt;/em&gt;). One particularly interesting trimming function is &lt;em&gt;robustness&lt;/em&gt; which is an extension of the concept of persistence to each points of the arcs ( an improved version of the &lt;em&gt;separatrix persistence&lt;/em&gt; described in &lt;em&gt;Weinkauf, T. and Gunther, D., 2009&lt;/em&gt; ). Indeed, &lt;em&gt;robustness&lt;/em&gt; can be considered a measure of how contrasted the filament is with respect to its &lt;em&gt;local&lt;/em&gt; background, and it is therefore a good way to select only highly contrasted subsets of the filamentary structures. Note that robustness, like persistence, is defined as a difference or ratio between the value of the field in two points, so the robustness threshold has the same order of magnitude as the persistence threshold.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;toRaDecZ&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-toRaDecZ&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Converts the coordinates to Ra (Right ascension), Dec (Declination) and Z (redshift). This is useful when the input file was computed from a particle distribution whose coordinates where given in the same system (such as a galaxy catalog for instance, see catalog format in &lt;a href=&quot;http://localhost/dotclear/index.php?post/field-formats&quot;&gt;field format&lt;/a&gt;).&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;toRaDecDist&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-toRaDecDist&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Converts the coordinates to Ra (Right ascension), Dec (Declination) and Dist (Distance). This is useful when the input file was computed from a particle distribution whose coordinates where given in the same system (such as a galaxy catalog for instance, see catalog format in &lt;a href=&quot;http://localhost/dotclear/index.php?post/field-formats&quot;&gt;field format&lt;/a&gt;).&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;rmBoundary&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-rmBoundary&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Remove the arcs and nodes that lay on the boundary or outside the domain of definition (e.g. nodes/arcs at infinity generated by boundary conditions).&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;rmOutside&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-rmOutside&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Remove the artificial arcs and nodes generated by the boundary conditions (e.g arcs/nodes at infinity) but keep boundaries (this is less restricitve than &lt;em&gt;-&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#rmBoundary&quot;&gt;rmBoundary&lt;/a&gt;&lt;/em&gt;).&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;toFITS&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-toFITS [&amp;lt;Xres&amp;gt;] [&amp;lt;Yres&amp;gt;] [&amp;lt;Zres&amp;gt;]&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Samples the skeleton to a FITS image with specified resolution. If the resolution is omitted (parameters &lt;em&gt;Xres&lt;/em&gt;, &lt;em&gt;Yres&lt;/em&gt; and/or &lt;em&gt;Zres&lt;/em&gt;), then the output file will have the same resolution as the network that was used to compute the skeleton.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;info&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-info&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Prints information on the input file, such as the number of arcs, critical points and the name and type of additional fields.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;addField&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-addField &amp;lt;filename&amp;gt; &amp;lt;field_name&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Tags each segment (i.e. pieces of arcs) and node (i.e. critical point) of the skeleton with the interpolated value of a field. The parameter &lt;em&gt;filename&lt;/em&gt; is the name of a readable &lt;a href=&quot;http://localhost/dotclear/index.php?post/field-formats&quot;&gt;regular grid field format&lt;/a&gt; file containing the grid to be interpolated, and &lt;em&gt;field_name&lt;/em&gt; is the name of the additional field in the output file.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;to&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-to &amp;lt;format&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Outputs a file in the selected writable &lt;a href=&quot;http://localhost/dotclear/index.php?post/skeleton-formats&quot;&gt;skeleton format&lt;/a&gt;. A list of possible parameter values can be obtained by running &lt;em&gt;skelconv&lt;/em&gt; without any argument. Note that skeleton files can be converted to &lt;a href=&quot;http://localhost/dotclear/index.php?post/networks&quot;&gt;NDnet format&lt;/a&gt; unstructured networks readable by &lt;a href=&quot;http://localhost/dotclear/index.php?post/netconv&quot;&gt;netconv&lt;/a&gt;.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;</content>
    
    

    
      </entry>
    
  <entry>
    <title>netconv</title>
    <link href="http://localhost/dotclear/index.php?post/netconv" rel="alternate" type="text/html"
    title="netconv" />
    <id>urn:md5:20ce034fb0fb466a86b6cac343e6a793</id>
    <published>2240-05-09T06:49:00+01:00</published>
    <updated>2012-10-02T16:20:01+02:00</updated>
    <author><name>thierry sousbie</name></author>
        <dc:subject>Manual</dc:subject>
            
    <content type="html">    &lt;p&gt;&lt;strong&gt;netconv&lt;/strong&gt;  is used to post-treat, convert and display information about  &lt;a href=&quot;http://localhost/dotclear/index.php?post/networks&quot;&gt;unstructured network format&lt;/a&gt; files. Run &lt;em&gt;netconv&lt;/em&gt; without argument to see a list of available input and output formats.&lt;/p&gt;


&lt;p&gt;&lt;br /&gt;
&lt;br /&gt;&lt;/p&gt;


&lt;p&gt;&lt;ins&gt;Usage:&lt;/ins&gt;&lt;/p&gt;

&lt;pre&gt; &lt;strong&gt;netconv&lt;/strong&gt; &amp;lt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/netconv#filename&quot;&gt;filename&lt;/a&gt;&amp;gt; 
         [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/netconv#outName&quot;&gt;outName&lt;/a&gt; &amp;lt;output filename&amp;gt;] [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/netconv#outDir&quot;&gt;outDir&lt;/a&gt; &amp;lt;dir&amp;gt;]
         [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/netconv#addField&quot;&gt;addField&lt;/a&gt; &amp;lt;filename&amp;gt; &amp;lt;field name&amp;gt;] [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/netconv#smooth&quot;&gt;smooth&lt;/a&gt; &amp;lt;Ntimes=0&amp;gt;] 
         [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/netconv#smoothData&quot;&gt;smoothData&lt;/a&gt; &amp;lt;vertex field name&amp;gt; &amp;lt;Ntimes&amp;gt;] 
         [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/netconv#noTags&quot;&gt;noTags&lt;/a&gt;] [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/netconv#toRaDecZ&quot;&gt;toRaDecZ&lt;/a&gt;] [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/netconv#toRaDecDist&quot;&gt;toRaDecDist&lt;/a&gt;] 
         [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/netconv#info&quot;&gt;info&lt;/a&gt;] [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/netconv#to&quot;&gt;to&lt;/a&gt; &amp;lt;format&amp;gt;]&lt;/pre&gt;

&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;filename&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;&amp;lt;filename&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;The name of a file containing an unstructured network (for instance, persistence pairs or manifolds as output by &lt;a href=&quot;http://localhost/dotclear/index.php?post/mse&quot;&gt;mse&lt;/a&gt;) in a readable &lt;a href=&quot;http://localhost/dotclear/index.php?post/networks&quot;&gt;network file format&lt;/a&gt;.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;outName&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-outName &amp;lt;fname&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Specifies the base name of the output file (extensions are added to this base name depending on the output type).&lt;br /&gt;&lt;ins&gt;Default value&lt;/ins&gt;: the name of the input file.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;outDir&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-outDir &amp;lt;dir&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Specifies the output directory. &lt;br /&gt;&lt;ins&gt;Default value&lt;/ins&gt;: the current working directory.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;addField&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-addField &amp;lt;filename&amp;gt; &amp;lt;field_name&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Tags each vertex of the network with the interpolated value of a field. The parameter &lt;em&gt;filename&lt;/em&gt; is the name of a readable &lt;a href=&quot;http://localhost/dotclear/index.php?post/field-formats&quot;&gt;regular grid field format&lt;/a&gt; file containing the grid to be interpolated, and &lt;em&gt;field_name&lt;/em&gt; is the name of the additional field in the output file.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;smooth&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-smooth &amp;lt;Ntimes=0&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Smooth the network N times. Smoothing is achieved by averaging the position of each vertex with that of its direct neighbors (i.e. those that belong to the boundary of at least one common simplex).  In practice, smoothing &lt;em&gt;N&lt;/em&gt; times makes the network smooth over the local size of '~N' network cells. Smoothing is first achieved only on vertices sharing at least a common 3-simplex, then a 2-simplex, ... So for instance when smoothing manifolds such as output by &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse&quot;&gt;mse&lt;/a&gt; -dumpManifolds JE0a&lt;/em&gt;, walls will be smoothed independently of filaments and critical points (try it yourself to understand ...).&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;smoothData&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-smoothData &amp;lt;vertex_field_name&amp;gt; &amp;lt;Ntimes&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;smooth the &lt;em&gt;vertex_field_name&lt;/em&gt; data field associated with vertices by averaging its value with that of its direct neighbors in the network &lt;em&gt;Ntimes&lt;/em&gt; times.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;noTags&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-noTags&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Prevents &lt;em&gt;netconv&lt;/em&gt; from adding trailing extensions to the output filename. Note that the last extension correponding to the file format is still added.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;toRaDecZ&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-toRaDecZ&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Converts the coordinates to Ra (Right ascension), Dec (Declination) and Z (redshift). This is useful when the input file was computed from a particle distribution whose coordinates where given in the same system (such as a galaxy catalog for instance, see catalog format in &lt;a href=&quot;http://localhost/dotclear/index.php?post/field-formats&quot;&gt;field format&lt;/a&gt;).&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;toRaDecDist&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-toRaDecDist&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Converts the coordinates to Ra (Right ascension), Dec (Declination) and Dist (Distance). This is useful when the input file was computed from a particle distribution whose coordinates where given in the same system (such as a galaxy catalog for instance, see catalog format in &lt;a href=&quot;http://localhost/dotclear/index.php?post/field-formats&quot;&gt;field format&lt;/a&gt;).&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;info&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-info&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Prints information on the input file, such as the number of each type of cell and the name and type of additional fields.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;to&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-to &amp;lt;format&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Outputs a file in the selected writable &lt;a href=&quot;http://localhost/dotclear/index.php?post/networks&quot;&gt;unstructured network format&lt;/a&gt;. A list of possible parameter values can be obtained by running &lt;em&gt;netconv&lt;/em&gt; without any argument.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;</content>
    
    

    
      </entry>
    
  <entry>
    <title>pdview</title>
    <link href="http://localhost/dotclear/index.php?post/pdview" rel="alternate" type="text/html"
    title="pdview" />
    <id>urn:md5:5c0c6c89724e22505ca34b99d696d65e</id>
    <published>2230-05-09T06:48:00+01:00</published>
    <updated>2012-06-03T16:41:57+02:00</updated>
    <author><name>thierry sousbie</name></author>
        <dc:subject>Manual</dc:subject>
            
    <content type="html">    &lt;p&gt;&lt;strong&gt;pdview&lt;/strong&gt; is used to interactively select the persistence threshold with the help of a persistence diagram. A persistence diagram is a 2D plot in which all the persistence pairs are represented by points with coordinates the value at the critical points in the pair (see &lt;a href=&quot;http://localhost/dotclear/index.php?post/Persistence-and-simplification&quot;&gt;this section&lt;/a&gt; for more info). In &lt;em&gt;pdview&lt;/em&gt;, we rather represent the value of the lower critical index on the x axis and the absolute difference in value between the two critical points on the y axis (i.e. the persistence of the pair). The diagram is computed from a &lt;a href=&quot;http://localhost/dotclear/index.php?post/networks&quot;&gt;network&lt;/a&gt; file encoding the persistence pairs as output by &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse&quot;&gt;mse&lt;/a&gt; filename -ppairs&lt;/em&gt;. Press the button labeled &lt;em&gt;?&lt;/em&gt; in &lt;strong&gt;pdview&lt;/strong&gt; for more information on how to use it.
&lt;br /&gt;
&lt;br /&gt;&lt;/p&gt;


&lt;p&gt;&lt;ins&gt;Usage:&lt;/ins&gt;&lt;/p&gt;

&lt;pre&gt; &lt;strong&gt;pdview&lt;/strong&gt; &amp;lt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/pdview#filename&quot;&gt;filename&lt;/a&gt;&amp;gt; [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/pdview#dtfe&quot;&gt;dtfe&lt;/a&gt;] [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/pdview#nsig&quot;&gt;nsig&lt;/a&gt; &amp;lt;val&amp;gt;] [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/pdview#cut&quot;&gt;cut&lt;/a&gt; &amp;lt;val&amp;gt;]&lt;/pre&gt;

&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;filename&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;&amp;lt;filename&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;The name of an &lt;a href=&quot;http://localhost/dotclear/index.php?post/networks&quot;&gt;unstructured network&lt;/a&gt; file containing the persistence pairs, as output by ''&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse&quot;&gt;mse&lt;/a&gt;.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;dtfe&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-dtfe&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Use this option if the persistence pairs were computed from a network estimated through DTFE (i.e. the output of &lt;a href=&quot;http://localhost/dotclear/index.php?post/Usage&quot;&gt;delaunay_nD&lt;/a&gt;).&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;nsig&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-nsig &amp;lt;val&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Sets a default persistence threshold ratio expressed in &lt;em&gt;N-sigmas&lt;/em&gt; (see option &lt;em&gt;-nsig&lt;/em&gt; in &lt;a href=&quot;http://localhost/dotclear/index.php?post/mse&quot;&gt;mse&lt;/a&gt;).&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;cut&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-cut &amp;lt;val&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Sets a default persistence threshold (see option &lt;em&gt;-cut&lt;/em&gt; in &lt;a href=&quot;http://localhost/dotclear/index.php?post/mse&quot;&gt;mse&lt;/a&gt;).&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;</content>
    
    

    
      </entry>
    
  <entry>
    <title>delaunay_nD</title>
    <link href="http://localhost/dotclear/index.php?post/Usage" rel="alternate" type="text/html"
    title="delaunay_nD" />
    <id>urn:md5:78515d55b2e404456778c4f305476e5f</id>
    <published>2220-05-09T05:22:00+01:00</published>
    <updated>2012-10-02T16:50:10+02:00</updated>
    <author><name>thierry sousbie</name></author>
        <dc:subject>Manual</dc:subject>
            
    <content type="html">    &lt;p&gt;&lt;strong&gt;Delaunay_2D&lt;/strong&gt; and &lt;strong&gt;Delaunay_3D&lt;/strong&gt; are used to compute the Delaunay tessellation of a discrete particle distribution (such as a N-Body simulation or a catalog of object coordinates) in 2D and 3D respectively and possibly in parallel. The output is an &lt;a href=&quot;http://localhost/dotclear/index.php?post/networks&quot;&gt;unstructured network&lt;/a&gt; in &lt;em&gt;NDnet&lt;/em&gt; format, with density computed for each vertex using DTFE (i.e. density at a given vertex is proportional to the total volume of the surrounding cells). In order to compute the correct topology and density close to the boundaries, the distribution is extrapolated over a band outside the domain of definition (see options &lt;em&gt;-&lt;a href=&quot;http://localhost/dotclear/index.php?post/Usage#margin&quot;&gt;margin&lt;/a&gt;&lt;/em&gt;, &lt;em&gt;-&lt;a href=&quot;http://localhost/dotclear/index.php?post/Usage#btype&quot;&gt;btype&lt;/a&gt;&lt;/em&gt; and &lt;em&gt;-&lt;a href=&quot;http://localhost/dotclear/index.php?post/Usage#periodic&quot;&gt;periodic&lt;/a&gt;&lt;/em&gt;). The output network can be used as input file for &lt;a href=&quot;http://localhost/dotclear/index.php?post/mse&quot;&gt;mse&lt;/a&gt; to compute the Morse-Smale complex of the initial discrete sample (the index of the vertices in the output network correspond to the index of the sampling particles in the input file).
&lt;br /&gt;
&lt;br /&gt;&lt;/p&gt;


&lt;p&gt;&lt;ins&gt;Usage:&lt;/ins&gt;&lt;/p&gt;

&lt;pre&gt; &lt;strong&gt;delaunay_nD&lt;/strong&gt; &amp;lt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/Usage#fname&quot;&gt;fname&lt;/a&gt;&amp;gt; [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/Usage#outName&quot;&gt;outName&lt;/a&gt; &amp;lt;fname&amp;gt;] [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/Usage#outDir&quot;&gt;outDir&lt;/a&gt; &amp;lt;dir&amp;gt;]
             [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/Usage#subbox&quot;&gt;subbox&lt;/a&gt; &amp;lt;x0&amp;gt; &amp;lt;y0&amp;gt; [&amp;lt;z0&amp;gt;] &amp;lt;dx&amp;gt; &amp;lt;dy&amp;gt; [&amp;lt;dz&amp;gt;]]
             [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/Usage#periodic&quot;&gt;periodic&lt;/a&gt;] [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/Usage#minimal&quot;&gt;minimal&lt;/a&gt;]
             [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/Usage#blocks&quot;&gt;blocks&lt;/a&gt; &amp;lt;NChunks&amp;gt; &amp;lt;NThreads&amp;gt;]
             [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/Usage#margin&quot;&gt;margin&lt;/a&gt; &amp;lt;M&amp;gt;] [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/Usage#btype&quot;&gt;btype&lt;/a&gt; &amp;lt;t=mirror&amp;gt;]
             [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/Usage#mask&quot;&gt;mask&lt;/a&gt; &amp;lt;fname.ND&amp;gt;]
             [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/Usage#smooth&quot;&gt;smooth&lt;/a&gt; &amp;lt;N=0&amp;gt;]
             [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/Usage#angmask&quot;&gt;angmask&lt;/a&gt; &amp;lt;fname.fits&amp;gt; [&amp;lt;maximum angular size (degrees) = 5.00&amp;gt;]]
             [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/Usage#radialDensity&quot;&gt;radialDensity&lt;/a&gt; &amp;lt;A&amp;gt; &amp;lt;Dr&amp;gt; &amp;lt;B&amp;gt;]
             [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/Usage#subSample&quot;&gt;subSample&lt;/a&gt; &amp;lt; 0&amp;lt;s&amp;lt;1 &amp;gt;]&lt;/pre&gt;

&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;fname&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;&amp;lt;fname&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;The name of a file containing the discrete particle coordinates in a &lt;a href=&quot;http://localhost/dotclear/index.php?post/field-formats&quot;&gt;field&lt;/a&gt; format.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;outName&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-outName &amp;lt;fname&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Specifies the base name of the file in output (extensions are added to this base name depending on the output type).&lt;br /&gt;&lt;ins&gt;Default value&lt;/ins&gt;: the name of the input  file.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;outDir&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-outDir &amp;lt;dir&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Specifies the output directory. &lt;br /&gt;&lt;ins&gt;Default value&lt;/ins&gt;: the current working directory.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;subbox&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-subbox &amp;lt;x0&amp;gt; &amp;lt;y0&amp;gt; [&amp;lt;z0&amp;gt;] &amp;lt;dx&amp;gt; &amp;lt;dy&amp;gt; [&amp;lt;dz&amp;gt;]&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Restricts the computation to a box of size [dx dy dz] and with origin [x0,y0,z0]. The box may be larger than the actual bounding box or cross its boundaries. The distribution outside the box will depend on the boundary type, as set with option &lt;em&gt;-&lt;a href=&quot;http://localhost/dotclear/index.php?post/Usage#btype&quot;&gt;btype&lt;/a&gt;&lt;/em&gt;. Note that the distribution within the margin of the subbox (see option &lt;em&gt;-&lt;a href=&quot;http://localhost/dotclear/index.php?post/Usage#margin&quot;&gt;margin&lt;/a&gt;&lt;/em&gt;) is still used to compute the topology close to the boundaries, so this option may be used to cut a large distribution into several small sub-boxes with matching boundary distribution.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;periodic&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-periodic&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Use this option to force periodic boundary conditions for the output network. When this option is used, &lt;em&gt;delaunay_nD&lt;/em&gt; will try to reconnect the simplices crossing the bounding box so that the space is compactified to a torus. The operation may fail if the margin size is too small, as cells crossing opposite boundaries may not match, so always check the output of the program as an estimation of the correct margin size will be indicated in that case (see also option &lt;em&gt;-&lt;a href=&quot;http://localhost/dotclear/index.php?post/Usage#margin&quot;&gt;margin&lt;/a&gt;&lt;/em&gt;). Note that this option wil set default boundary type to &lt;em&gt;-&lt;a href=&quot;http://localhost/dotclear/index.php?post/Usage#btype&quot;&gt;btype&lt;/a&gt; periodic&lt;/em&gt;.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;minimal&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-minimal&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;when using this option, only the minimal necessary information to define the tesselation is stored in the output file (i.e. In dimensions D, the vertex coordinates and simplices of dimension D, and the additional data associated to them). Note that it is preferable NOT to use this option when the ouput file is to be fed to &lt;a href=&quot;http://localhost/dotclear/index.php?post/mse&quot;&gt;mse&lt;/a&gt;, as it will force the program to recompute the topology of intermediate simplices using a slower algorithm every time it is run.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;blocks&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-blocks &amp;lt;NChunks&amp;gt; &amp;lt;NThreads&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;instead of computing the delaunay tesselation of the full distribution, divide it into &lt;em&gt;NChunks&lt;/em&gt; overlapping sub blocks and process them &lt;em&gt;NThreads&lt;/em&gt; at a time. The subblocks are then automatically reassembled into the full delaunay tesselation. This option can either be used to increase speed by parallelizing the process (for high values of &lt;em&gt;NThreads&lt;/em&gt;) or decrease the memory consumption (when &lt;em&gt;NChunks&lt;/em&gt;&amp;gt;&amp;gt;&lt;em&gt;NThreads&lt;/em&gt;). Note that it is incompatible with options &lt;em&gt;-&lt;a href=&quot;http://localhost/dotclear/index.php?post/Usage#mask&quot;&gt;mask&lt;/a&gt;&lt;/em&gt;, &lt;em&gt;-&lt;a href=&quot;http://localhost/dotclear/index.php?post/Usage#btype&quot;&gt;btype&lt;/a&gt; smooth&lt;/em&gt;, &lt;em&gt;-&lt;a href=&quot;http://localhost/dotclear/index.php?post/Usage#angmask&quot;&gt;angmask&lt;/a&gt;&lt;/em&gt; and &lt;em&gt;-&lt;a href=&quot;http://localhost/dotclear/index.php?post/Usage#radialDensity&quot;&gt;radialDensity&lt;/a&gt;&lt;/em&gt;.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;margin&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-margin &amp;lt;M&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Suggests a startup trial size for the additional margin used to compute the topology and density close to the boundary. The margin size is expressed as a fraction of the bounding box (or subbox if option &lt;em&gt;-&lt;a href=&quot;http://localhost/dotclear/index.php?post/Usage#subbox&quot;&gt;subbox&lt;/a&gt;&lt;/em&gt; is used). Note the the program can probably make a better guess than you so it is not recommanded to set it by hand ...&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;btype&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-btype &amp;lt;t=mirror&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;This option sets the type of the boundaries. This option is used to set how the distribution should be extrapolated  outside the bounding box (an estimation of the distribution outside the bounding box is needed to correctly estimate the topology and density of the distribution close to its boundaries). Possible boundary types are:
&lt;ul&gt;
&lt;li&gt;- &lt;em&gt;mirror&lt;/em&gt; : the distribution outside the bounding box is a mirrored copy of that inside.&lt;/li&gt;
&lt;li&gt;- &lt;em&gt;periodic&lt;/em&gt; : use periodic boundary condition (i.e. space is paved with copies of the distribution in the bounding box). Note that this option does *NOT* enforce periodic boundary conditions as it does not tell &lt;em&gt;delaunay_nD&lt;/em&gt; to reconnect the Delaunay cells that cross the bounding box (this is achieved with &lt;em&gt;-&lt;a href=&quot;http://localhost/dotclear/index.php?post/Usage#periodic&quot;&gt;periodic&lt;/a&gt;&lt;/em&gt;).&lt;/li&gt;
&lt;li&gt;- &lt;em&gt;smooth&lt;/em&gt; : a surface of guard particles is added outside the bounding box and new particles are added by interpolating the estimated density computed on the boundary of the distribution. This boundary type is useful when the actual boundaries of the sample are complex (i.e. not a cube), such as for a 3D galaxy catalog limited to a portion of the sky.&lt;/li&gt;
&lt;li&gt;- &lt;em&gt;void&lt;/em&gt; : the distribution is supposed to be void outside the bounding box.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;mask&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-mask &amp;lt;fname.ND&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Specifies a mask. The file must be a 1D array of values of size the number of vertices (or pixels) in the network, in a readable &lt;a href=&quot;http://localhost/dotclear/index.php?post/field-formats&quot;&gt;grid format&lt;/a&gt;. A value of &lt;em&gt;0&lt;/em&gt; corresponds to a non-masked particle while any other value masks the particle. Note that masked particles are still used to compute the Delaunay tessellation and estimate the density but they will change the topology of the actual distribution.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;smooth&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-smooth &amp;lt;N=0&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Smooth the estimated DTFE density by averaging N times its value at each vertex with that at the neighboring vertices (i.e two vertices are considered neighbors if they both belong to the boundary of at list one given cell).&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;angmask&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-angmask &amp;lt;fname.fits&amp;gt; [&amp;lt;maximum angular size (degrees) = 5.00&amp;gt;]&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Not documented ... you should not use it anyway :)&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;radialDensity&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-radialDensity &amp;lt;A&amp;gt; &amp;lt;Dr&amp;gt; &amp;lt;B&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Not documented ...&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;subSample&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-subSample &amp;lt; 0&amp;lt;s&amp;lt;1 &amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Computes the Delaunay tessellation over a subsets of the actual distribution containing only a randomly selected fraction &lt;em&gt;0&amp;lt;s&amp;lt;1&lt;/em&gt; of the particles.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;</content>
    
    

    
      </entry>
    
  <entry>
    <title>mse</title>
    <link href="http://localhost/dotclear/index.php?post/mse" rel="alternate" type="text/html"
    title="mse" />
    <id>urn:md5:eb9ed29b1b72c32779d712599fe5fb0d</id>
    <published>2210-05-09T05:21:00+01:00</published>
    <updated>2012-11-15T14:10:20+01:00</updated>
    <author><name>thierry sousbie</name></author>
        <dc:subject>Manual</dc:subject>
            
    <content type="html">    &lt;p&gt;&lt;strong&gt;mse&lt;/strong&gt; is the main program in DisPerSE, it is used to compute Morse-smale complexes and extract filamentary structures (i.e. arcs), voids and walls (ascending/descending manifolds), critical points and persistence pairs. By default, when no option is specified, the Morse-smale complex is computed and stored as a &lt;em&gt;filename.MSC&lt;/em&gt; backup file in a format internal to mse, and no other file is produced. This file may later be loaded by mse (option &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#loadMSC&quot;&gt;-loadMSC&lt;/a&gt;&lt;/em&gt;) to skip the computation of the Morse-smale complex (this is useful when you want to produce different outputs from the same source).
&lt;br /&gt;
&lt;br /&gt;&lt;/p&gt;


&lt;p&gt;&lt;ins&gt;Usage:&lt;/ins&gt;&lt;/p&gt;

&lt;pre&gt; &lt;strong&gt;mse&lt;/strong&gt; &amp;lt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#networkfilename&quot;&gt;network filename&lt;/a&gt;&amp;gt; [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#field&quot;&gt;field&lt;/a&gt; &amp;lt;fname&amp;gt;]
     [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#outname&quot;&gt;outName&lt;/a&gt; &amp;lt;fname&amp;gt;] [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#notags&quot;&gt;noTags&lt;/a&gt;] [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#outdir&quot;&gt;outDir&lt;/a&gt; &amp;lt;dir&amp;gt;] 
     [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#periodicity&quot;&gt;periodicity&lt;/a&gt; &amp;lt;val&amp;gt;]
     [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#mask&quot;&gt;mask&lt;/a&gt; &amp;lt;fname&amp;gt;[~]]
     [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#nthreads&quot;&gt;nthreads&lt;/a&gt; &amp;lt;N=$OMP_NUM_THREADS&amp;gt;]
 
     [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#nsig&quot;&gt;nsig&lt;/a&gt; &amp;lt;n1, n2, ...&amp;gt;] [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#cut&quot;&gt;cut&lt;/a&gt; &amp;lt;l1, l2, ...&amp;gt;]
     [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#interactive&quot;&gt;interactive&lt;/a&gt; [&amp;lt;full/path/to/pdview&amp;gt;]]
     [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#forceloops&quot;&gt;forceLoops&lt;/a&gt;]
     [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#robustness&quot;&gt;robustness&lt;/a&gt;] [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#robustness&quot;&gt;no_robustness&lt;/a&gt;]
 
     [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#manifolds&quot;&gt;manifolds&lt;/a&gt;] [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#interarcsgeom&quot;&gt;interArcsGeom&lt;/a&gt;] [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#no_arcsgeom&quot;&gt;no_arcsGeom&lt;/a&gt;]
 
     [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#ppairs&quot;&gt;ppairs&lt;/a&gt;] [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#ppairs_ASCII&quot;&gt;ppairs_ASCII&lt;/a&gt;] 
     [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#upskl&quot;&gt;upSkl&lt;/a&gt;] [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#downskl&quot;&gt;downSkl&lt;/a&gt;] [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#interskl&quot;&gt;interSkl&lt;/a&gt;]
     [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#dumparcs&quot;&gt;dumpArcs&lt;/a&gt; &amp;lt;CUID&amp;gt;]
     [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#dumpmanifolds&quot;&gt;dumpManifolds&lt;/a&gt; [&amp;lt;JEP0123456789ad&amp;gt;]] 
 
     [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#compactify&quot;&gt;compactify&lt;/a&gt; &amp;lt;type=natural&amp;gt;]
     [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#vertexasminima&quot;&gt;vertexAsMinima&lt;/a&gt;] [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#descendingfiltration&quot;&gt;descendingFiltration&lt;/a&gt;] 
     [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#no_savemsc&quot;&gt;no_saveMSC&lt;/a&gt;] [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#loadMSC&quot;&gt;loadMSC&lt;/a&gt; &amp;lt;fname&amp;gt;]
 
     [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#no_gfilter&quot;&gt;no_gFilter&lt;/a&gt;]
     [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#diagram&quot;&gt;diagram&lt;/a&gt;] [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#smooth&quot;&gt;smooth&lt;/a&gt; &amp;lt;Ntimes=0&amp;gt;] 
     [-&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#debug&quot;&gt;debug&lt;/a&gt;]&lt;/pre&gt;

&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;networkfilename&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;&amp;lt;network filename&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;The cell complex defining the topology of the space and optionally the function discretely sampled over it. The file may be a Healpix &lt;em&gt;FITS&lt;/em&gt; file, a &lt;a href=&quot;http://localhost/dotclear/index.php?post/field-formats&quot;&gt;regular grid&lt;/a&gt; or an &lt;a href=&quot;http://localhost/dotclear/index.php?post/networks&quot;&gt;unstructured network&lt;/a&gt; readable by &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/netconv&quot;&gt;netconv&lt;/a&gt;&lt;/em&gt; or &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/fieldconv&quot;&gt;fieldconv&lt;/a&gt;&lt;/em&gt; (run &lt;em&gt;netconv&lt;/em&gt; or &lt;em&gt;fieldconv&lt;/em&gt; without argument for a list of supported file formats). The value of the function from which the Morse-smale complex will be computed should be given for each vertex (or pixel in case of a regular grid). In the case of an unstructured network, the function should be given in a field labeled &lt;em&gt;field_value&lt;/em&gt;, or the function may also be set using option &lt;em&gt;-&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#field&quot;&gt;field&lt;/a&gt;&lt;/em&gt;.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;field&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-field &amp;lt;fname&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Specifies the scalar function whose morse complex will be computed. The function value should be given for each vertex (or pixel for regular grids), so the file must be a 1D array of values of size the number of vertices (or pixels) in the network, in a readable &lt;a href=&quot;http://localhost/dotclear/index.php?post/field-formats&quot;&gt;grid format&lt;/a&gt;.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;outname&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-outName &amp;lt;fname&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Specifies the base name of the output file (extensions are added to this base name depending on the output type).&lt;br /&gt;&lt;ins&gt;Default value&lt;/ins&gt;: the name of the input network file.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;outdir&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-outDir &amp;lt;dir&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Specifies the output directory. &lt;br /&gt;&lt;ins&gt;Default value&lt;/ins&gt;: the current working directory.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;notags&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-noTags&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Prevents mse from adding trailing extensions to the output filename such as the persistence cut levels ... Note that the last extension correponding to the file format is still added.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;periodicity&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-periodicity &amp;lt;val&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Specifies the periodicity of the domain. This is only applicable to grid types networks, which are considered by default to have periodic boundaries. The parameter &lt;em&gt;&amp;lt;val&amp;gt;&lt;/em&gt; is a serie of &lt;em&gt;1&lt;/em&gt;s and &lt;em&gt;0&lt;/em&gt;s enabling/disabling periodic boundary conditions along the corresponding direction.&lt;br /&gt;&lt;ins&gt;Exemple&lt;/ins&gt;: &lt;em&gt;-periodicity 0&lt;/em&gt; sets non periodic boundaries conditions (PBC) while &lt;em&gt;-periodicity 101&lt;/em&gt; sets PBC along dims &lt;em&gt;0&lt;/em&gt; and &lt;em&gt;2&lt;/em&gt; (&lt;em&gt;x&lt;/em&gt; and &lt;em&gt;z&lt;/em&gt;) but not along &lt;em&gt;y&lt;/em&gt; axis.&lt;br /&gt;&lt;ins&gt;Default value&lt;/ins&gt;: by default, boundary conditions are fully periodic.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;mask&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-mask &amp;lt;fname&amp;gt; [~]&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Specifies a mask. The file must be a 1D array of values of size the number of vertices (or pixels) in the network, in a readable &lt;a href=&quot;http://localhost/dotclear/index.php?post/field-formats&quot;&gt;grid format&lt;/a&gt;. By default, a value of &lt;em&gt;0&lt;/em&gt; corresponds to a visible vertex/pixel while any other value masks the vertex/pixel. Adding a trailing &lt;em&gt;~&lt;/em&gt; to the filename (without space) reverses this behavior, a value of 0 masking the corresponding pixels/vertices.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;nthreads&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-nthreads &amp;lt;N=$OMP_NUM_THREADS&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Specifies the number of threads.&lt;br /&gt;&lt;ins&gt;Default value&lt;/ins&gt;: the number of threads is set by the environment variable &lt;em&gt;$OMP_NUM_THREADS&lt;/em&gt;  which is usually set to the total number of cores available by &lt;em&gt;openMP&lt;/em&gt;.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;nsig&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-nsig &amp;lt;n1, n2, ...&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Specifies the persistence ratio threshold in terms of &quot;number of sigmas&quot;. Any persistence pair with a persistence ratio (i.e. the ratio of the values of the points in the pair) that has a probability less than &quot;n-sigmas&quot; to appear in a random field will be cancelled. This may only be used for discretely sampled density fields (such as N-body simulations or discrete objects catalogs) whose density was estimated through the DTFE of a delaunay tesselation. This option is typically used when the input network was produced using &lt;a href=&quot;http://localhost/dotclear/index.php?post/Usage&quot;&gt;delaunay_2D&lt;/a&gt; or &lt;a href=&quot;http://localhost/dotclear/index.php?post/Usage&quot;&gt;delaunay_3D&lt;/a&gt;, in any other case, use &lt;em&gt;-cut&lt;/em&gt; instead. See also option &lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#interactive&quot;&gt;-interactive&lt;/a&gt;&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;cut&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-cut &amp;lt;l1, l2, ...&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Specifies the persistence threshold. Any persistence pair with persistence (i.e. the difference of value between the two points in the pair) lower than the given threshold will be cancelled. Use &lt;em&gt;-nsig&lt;/em&gt; instead of this when the input network was produced with &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/Usage&quot;&gt;delaunay_2D&lt;/a&gt;&lt;/em&gt; or &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/Usage&quot;&gt;delaunay_3D&lt;/a&gt;&lt;/em&gt;. The cut value should be typically set to the estimated amplitude of the noise. See also option &lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#interactive&quot;&gt;-interactive&lt;/a&gt;.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;interactive&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-interactive [&amp;lt;full/path/to/pdview&amp;gt;]&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;This option allows the user to interactively select the persistence threshold on a persistence diagram (using &lt;a href=&quot;http://localhost/dotclear/index.php?post/pdview&quot;&gt;pdview&lt;/a&gt;). if option &lt;em&gt;-cut&lt;/em&gt; or &lt;em&gt;-nsig&lt;/em&gt; is also used, the specified values are become default thresholds in pdview. The path to &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/pdview&quot;&gt;pdview&lt;/a&gt;&lt;/em&gt; may be optionally indicated in case it cannot be found automatically.&lt;br /&gt;&lt;ins&gt;Note&lt;/ins&gt;: this option may only be used if pdview was compiled (requires &lt;em&gt;mathGL&lt;/em&gt; and &lt;em&gt;Qt4&lt;/em&gt; libraries).&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;forceloops&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-forceLoops&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Forces the simplification of non-cancellable persistence pairs (saddle-saddle pairs in 3D or more that are linked by at least 2 different arcs). When two critical points of critical index difference 1 are linked by 2 or more arcs, they may not be cancelled as this would result in a discrete gradient loop. This is not a problem in 2D as such pairs cannot form persistence pairs but in 3D, saddle-saddle persistence pairs may be linked by 2 or more arcs even though their persistence is low. By default those pairs are skipped in order to preserve the properties of the Morse-smale complex but as a result few non persistent features may remain (such as spurious filaments). Fortunately, there are usually none or only very few of those pairs, and their number is shown in the output of mse, in the &lt;em&gt;Simplifying complex&lt;/em&gt; section. If you are only interested in identifying structures (as opposed to doing topology),  you should probably use '-forceLoops' to remove those spurious structures (such as small non significant filaments).&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;robustness&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-robustness / -no_robustness&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Enables or prevents the computation of &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/definitions#robustness&quot;&gt;robustness&lt;/a&gt;&lt;/em&gt; and &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/definitions#robustness&quot;&gt;robustness ratio&lt;/a&gt;&lt;/em&gt;. By default, &lt;em&gt;robustness&lt;/em&gt; is not computed as it can be costly for very large data sets. When enabled, a robustness value is tagged for each segments and node of the output &lt;a href=&quot;http://localhost/dotclear/index.php?post/skeleton-formats&quot;&gt;skeleton&lt;/a&gt; files. See also options &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#trim&quot;&gt;-trimBelow&lt;/a&gt;&lt;/em&gt; of &lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv&quot;&gt;skelconv&lt;/a&gt; and the &lt;a href=&quot;http://localhost/dotclear/index.php?post/regular-grid-%3A-filaments#robustness&quot;&gt;tutorial section&lt;/a&gt; for applications.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;manifolds&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-manifolds&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Forces the computation and storage of all &lt;a href=&quot;http://localhost/dotclear/index.php?post/definitions&quot;&gt;ascending and descending manifolds&lt;/a&gt; (i.e. walls, voids, ...). By default, mse only stores manifolds geometry if required (for instance, when the option &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#dumpmanifolds&quot;&gt;-dumpManifolds&lt;/a&gt;&lt;/em&gt; is used), and the resulting backup MSC file will therefore not contain this information and may not be used later to compute manifolds.&lt;br /&gt;&lt;ins&gt;Example&lt;/ins&gt;: running the command 'mse filename' will produce a the backup file 'filename.MSC' that stores the Morse-smale complex information. A later call to &lt;em&gt;mse -loadMSC filename.MSC -upSkl&lt;/em&gt; would skip the computation of the MS-complex and succeed in producing the skeleton of the filamentary structures (arcs) as arcs geometry is computed by default. Unfortunately, &lt;em&gt;mse -loadMSC filename.MSC -dumpManifolds J0a&lt;/em&gt; would fail as the information needed to compute the ascending 3-manifolds (i.e. manifolds originating from a minimum, which trace the voids) is not avaialble by default. Running &lt;em&gt;mse filename -manifolds&lt;/em&gt; in the beginning would solve this problem.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;interarcsgeom&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-interArcsGeom&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;This is similar to option &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#manifolds&quot;&gt;-manifolds&lt;/a&gt;&lt;/em&gt;, but for the arcs linking different types of saddle-points together (in 3D or more). By default, unless option &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#interskl&quot;&gt;-interSkl&lt;/a&gt;&lt;/em&gt; is specified, only the geometry of arcs linking extrema to saddle points is computed, so the backup &lt;em&gt;.MSC&lt;/em&gt; file may not be used later to retrieve other types of arcs geometry. Specifying &lt;em&gt;-interArcsGeom&lt;/em&gt; forces the computation and storage of all types of arcs in the &lt;em&gt;.MSC&lt;/em&gt; backup file.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;no_arcsgeom&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-no_arcsGeom&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;By default, the geometry of arcs linking extrema and saddle point is always computed even when not directly needed, so that it is stored in the backup &lt;em&gt;.MSC&lt;/em&gt; file for later use. Specifying &lt;em&gt;-no_arcsGeom&lt;/em&gt; may be used to lower memory usage when only the critical points and/or their persistence pairings are needed.&lt;br /&gt;&lt;ins&gt;Exemple&lt;/ins&gt;: &lt;em&gt;mse filename -ppairs -no_arcsGeom&lt;/em&gt;  will compute the critical points and persistence pairs using less memory than &lt;em&gt;mse filename -ppairs&lt;/em&gt;, but the resulting &lt;em&gt;.MSC&lt;/em&gt; file may not be later used to retrieve arcs geometry.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;ppairs&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-ppairs&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Dumps the persistence pairs as a &lt;em&gt;NDnet&lt;/em&gt; &lt;a href=&quot;http://localhost/dotclear/index.php?post/networks&quot;&gt;network type file&lt;/a&gt;. The resulting file is a 1D network, where the critical points are the vertices and 1-cells (segments) represent pairs. Additional information such as the type, persistence, cell in the initial complex, ... of each critical points and pairs is also stored as additional information. Try running &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/netconv&quot;&gt;netconv&lt;/a&gt; filename.ppairs.NDnet -info&lt;/em&gt; for a list of available additional data (see also &lt;em&gt;additional data&lt;/em&gt; in the &lt;a href=&quot;http://localhost/dotclear/index.php?post/networks&quot;&gt;network file format&lt;/a&gt; section).&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;ppairs_ascii&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-ppairs_ASCII&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Same as &lt;em&gt;-&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#ppairs&quot;&gt;ppairs&lt;/a&gt;&lt;/em&gt; but pairs are dumped in a easily readable ASCII format. This option is deprecated and should not be used, use &lt;em&gt;-ppairs&lt;/em&gt; option instead and then run &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/netconv&quot;&gt;netconv&lt;/a&gt; filename.ppairs.NDnet -to NDnet_ascii&lt;/em&gt;.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;upskl&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-upSkl&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Dumps the &quot;up&quot; skeleton (i.e. arcs linking maxima to saddle-points, which trace the filamentary structures) as a &lt;em&gt;NDskl&lt;/em&gt; type &lt;a href=&quot;http://localhost/dotclear/index.php?post/skeleton-formats&quot;&gt;skeleton file&lt;/a&gt;. This command is an alias of &lt;em&gt;-&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#dumparcs&quot;&gt;dumpArcs&lt;/a&gt; U&lt;/em&gt;.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;downskl&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-downSkl&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Dumps the &quot;down&quot; skeleton (i.e. arcs linking minima to saddle-points, which trace the anti-filaments, or void filaments) as a &lt;em&gt;NDskl&lt;/em&gt; type &lt;a href=&quot;http://localhost/dotclear/index.php?post/skeleton-formats&quot;&gt;skeleton file&lt;/a&gt;. This command is an alias of &lt;em&gt;-&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#dumparcs&quot;&gt;dumpArcs&lt;/a&gt; D&lt;/em&gt;.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;interskl&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-interSkl&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Dumps the &quot;inter&quot; skeleton (i.e. arcs linking different types of saddle-points together) as a &lt;em&gt;NDskl&lt;/em&gt; type &lt;a href=&quot;http://localhost/dotclear/index.php?post/skeleton-formats&quot;&gt;skeleton file&lt;/a&gt;. This will only work in 3D or more, and be careful that those type of arcs may have a very complex structure. This command is an alias of '-&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#dumparcs&quot;&gt;dumpArcs&lt;/a&gt; I''.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;dumparcs&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-dumpArcs &amp;lt;CUID&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Dumps the specified type of &lt;a href=&quot;http://localhost/dotclear/index.php?post/definitions&quot;&gt;arcs&lt;/a&gt; (i.e filaments) as a &lt;em&gt;NDskl&lt;/em&gt; type &lt;a href=&quot;http://localhost/dotclear/index.php?post/skeleton-formats&quot;&gt;skeleton file&lt;/a&gt;. Any combination of the letter &lt;em&gt;C&lt;/em&gt;, &lt;em&gt;U&lt;/em&gt;, &lt;em&gt;I&lt;/em&gt; and &lt;em&gt;D&lt;/em&gt; may be used as parameter to indicate which type of arc should be saved:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;em&gt;U&lt;/em&gt;(p): arcs leading to maxima (see also &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#upskl&quot;&gt;-upSkl&lt;/a&gt;&lt;/em&gt;).&lt;/li&gt;
&lt;li&gt;-&lt;em&gt;D&lt;/em&gt;(own): arcs leading to minima (see also &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#downskl&quot;&gt;-downSkl&lt;/a&gt;&lt;/em&gt;).&lt;/li&gt;
&lt;li&gt;-&lt;em&gt;I&lt;/em&gt;(nter): other arcs linking saddle-points together (see also &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#interskl&quot;&gt;-interSkl&lt;/a&gt;&lt;/em&gt;).&lt;/li&gt;
&lt;li&gt;-&lt;em&gt;C&lt;/em&gt;(onnect): keeps at least the connectivity information for all types of arcs.&lt;br /&gt;&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;li&gt;&lt;ins&gt;Exemple&lt;/ins&gt;: &lt;em&gt;CUD&lt;/em&gt; dumps geometry of arcs for which one extremity is and extremum (a maximum or a minimum), and the connectivity of saddle points is also stored (i.e. one can retrieve how they are linked in the Morse-Smale complex but not the actual geometry of the corresponding arcs).&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;dumpmanifolds&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-dumpManifolds [&amp;lt;JEP0123456789ad&amp;gt;]&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Dumps &lt;a href=&quot;http://localhost/dotclear/index.php?post/definitions&quot;&gt;ascending and/or descending manifolds&lt;/a&gt; (i.e. walls, voids, ...) as a &lt;em&gt;NDnet&lt;/em&gt; type &lt;a href=&quot;http://localhost/dotclear/index.php?post/networks&quot;&gt;network file&lt;/a&gt;. The first part of the argument must be a combination of letters &lt;em&gt;J&lt;/em&gt;, &lt;em&gt;E&lt;/em&gt; and &lt;em&gt;P&lt;/em&gt; and the second, which indicates the type of manifold to dump as a digit indicating the critical index of the critical point it originates from followed by the letter &lt;em&gt;a&lt;/em&gt; and/or &lt;em&gt;d&lt;/em&gt; to select ascending and/or descending manifolds respectively:
&lt;ul&gt;
&lt;li&gt;- &lt;em&gt;J&lt;/em&gt;(oin): join all the the manifolds in one single file. By default, one file is created for each manifold.&lt;/li&gt;
&lt;li&gt;- &lt;em&gt;E&lt;/em&gt;(xtended): compute extended manifolds. For instance, when computing an ascending 2-manifold (i.e. a wallseparating two voids), the ascending 1-manifolds (i.e. filaments) on its boundary as well as the ascending 0-manifolds (i.e. critical points) on their boundaries are also dumped.&lt;/li&gt;
&lt;li&gt;- &lt;em&gt;P&lt;/em&gt;(reserve): do not merge infinitely close submanifolds. By default, the overlapping part of manifolds are merged (e.g. filaments may have bifurcation points but arcs only stop at critical points, so although two arcs continue above the bifurcation, only one path is stored). Using this option preserve the properties of the MS-complex but may result in very large files.&lt;/li&gt;
&lt;li&gt;- &lt;em&gt;D&lt;/em&gt;(ual): compute dual cells geometry when appropriate. Descending and ascending manifolds are the dual of one another, so they should be represented by cells of the complex and their dual respectively. When this option is used, appropriate cells are replaced by their dual (for instance, ascending 0-manifolds of a delaunay tesselation are represented by sets of voronoi 3-cells, the dual of the 0-cells (vertices) of the delaunay tesselation). The default behavior is that descending manifolds belong to the dual, so ascending manifolds are always fine but one may want to use option &lt;em&gt;D&lt;/em&gt; to compute dual cells or option &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#vertexasminima&quot;&gt;-vertexAsMinima&lt;/a&gt;&lt;/em&gt; to change this behavior in order to obtain correct descending manifolds.&lt;/li&gt;
&lt;li&gt;- &lt;em&gt;0123456789&lt;/em&gt;: specifies the critical index of the critical points from which the manifold originates (ascending 0,1,2 and 3-manifolds in 3D trace the voids, walls, filaments and maxima/sources respectively).&lt;/li&gt;
&lt;li&gt;- &lt;em&gt;a/d&lt;/em&gt;: compute the ascending and/or descending manifold.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;li&gt;&lt;ins&gt;Exemple 1&lt;/ins&gt;: in 3D, argument &lt;em&gt;JD0a2ad3d&lt;/em&gt; would dump to a single file the ascending manifolds of critical index 0 and 2, and the descending ones for critical index 2 and 3. Cells are replaced by their dual where appropriate (i.e. in 3D, the descending 2-manifolds would be 2-cells dual of segments, and descending 3 manifolds would be 3-cells dual to vertices).&lt;/li&gt;
&lt;li&gt;&lt;ins&gt;Exemple 2&lt;/ins&gt;: in 3D, argument &lt;em&gt;JE0a&lt;/em&gt; would dump the union of the tetrahedrons within the voids, the walls separating the voids, as well as the filaments on their boundaries and the maxima at their extremity, as sets of triangles, segments and vertices of the input complex respectively, to a single file. As option &lt;em&gt;P&lt;/em&gt; was not used, any overlapping region would be stored only once.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;compactify&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-compactify &amp;lt;type=natural&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;This option is used to specify how to treat boundaries of manifolds with boundary (this option does not affect manifolds without boundary). Available type are:
&lt;ul&gt;
&lt;li&gt;- 'natural' is the default value, and almost always the best choice.&lt;/li&gt;
&lt;li&gt;- 'torus' creates reflective (i.e periodic) boundaries.&lt;/li&gt;
&lt;li&gt;- 'sphere' links boundaries to a cell at infinity.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;vertexasminima&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-vertexAsMinima&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;As mse uses discrete Morse theory, each type of critical point corresponds to a particular cell type. By defaults, 0-cells (vertices) are maxima, 1-cells are saddle points, .... and d-cells are minima (d is the number of dimensions). Using &lt;em&gt;-vertexAsMinima&lt;/em&gt;, the association is reversed and 0-cells are associated to minima ... while d-cells are associated to maxima. This option can be useful to identify manifolds or arcs as lists of a particular type of cell.&lt;br /&gt;&lt;ins&gt;Example&lt;/ins&gt;: &lt;em&gt;mse input_filename -&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#dumpmanifolds&quot;&gt;dumpManifolds&lt;/a&gt; J0a -vertexAsMinima&lt;/em&gt; can be used to identify voids (ascending 0-manifolds) as sets of 0-cells (vertices). As a result, voids are identified as sets of pixels or  particles (if input_filename is a grid or the delaunay tesselation of an N-body simulation respectiveley), which is easy to relate to the pixels/particles of the input file. If the command &lt;em&gt;mse input_filename -dumpManifolds J0a&lt;/em&gt; had been issued instead, each void would have been described by a set of n-cells (in nD).&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;descendingfiltration&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-descendingFiltration&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;By default, &lt;em&gt;mse&lt;/em&gt; uses an ascending filtration to compute the discrete gradient. This option forces the program to use a descending filtration instead. Not that if mse works correctly (and it should hopefully be the case) the general properties of the Morse-smale complex are not affected, but simply critical points and  manifolds geometry will slightly differ.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;no_savemsc&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-no_saveMSC&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Use this if you do not want &lt;em&gt;mse&lt;/em&gt; to write a backup &lt;em&gt;.MSC&lt;/em&gt; file.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;loadMSC&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-loadMSC &amp;lt;fname&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Loads a given backup &lt;em&gt;.MSC&lt;/em&gt; file. This will basically skip the computation of the Morse-smale complex, therefore gaining a lot of time. By default, mse always write a backup &lt;em&gt;.MSC&lt;/em&gt; after the MS-complex is computed. See options &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#manifolds&quot;&gt;-manifolds&lt;/a&gt;&lt;/em&gt; and &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#interarcsgeom&quot;&gt;-interArcsGeom&lt;/a&gt;&lt;/em&gt; for information on the limitations.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;no_gfilter&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-no_gFilter&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Prevents the filtration of null-persistence pairs in the discrete gradient. There is no point in using this option except for debugging ...&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;diagram&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-diagram&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Draws a diagram of the connectivity of critical points to a &lt;em&gt;.ps&lt;/em&gt; file. This should only be used when the number of critical points is very small (i.e. less than 100 or so).&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;smooth&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-smooth &amp;lt;Ntimes=0&amp;gt;&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Deprecated, don't use it ... use &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv&quot;&gt;skelconv&lt;/a&gt; -smooth &amp;lt;Ntimes=0&amp;gt;&lt;/em&gt; or &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/netconv&quot;&gt;netconv&lt;/a&gt; -smooth &amp;lt;Ntimes=0&amp;gt;&lt;/em&gt; on the output files instead.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;debug&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;&lt;strong&gt;-debug&lt;/strong&gt;:&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Outputs some debug information and files.&lt;/li&gt;
&lt;/ul&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;</content>
    
    

    
      </entry>
  
</feed>