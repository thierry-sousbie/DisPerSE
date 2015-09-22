<?xml version="1.0" encoding="utf-8"?><feed xmlns="http://www.w3.org/2005/Atom"
  xmlns:dc="http://purl.org/dc/elements/1.1/"
  xmlns:wfw="http://wellformedweb.org/CommentAPI/"
  xml:lang="en">
  
  <title type="html">DisPerSE - persistent structures identification - Skeleton data</title>
  <subtitle type="html">Automatic identification of persistent structures in 2D or 3D.
keywords: Morse complex, topology, peak, void, source, wall, filament, cosmic web, cosmology, structure identification.</subtitle>
  <link href="http://localhost/dotclear/index.php?feed/category/Skeleton-I-O/atom" rel="self" type="application/atom+xml"/>
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
    <title>vtk skeleton formats</title>
    <link href="http://localhost/dotclear/index.php?post/vtk-formats" rel="alternate" type="text/html"
    title="vtk skeleton formats" />
    <id>urn:md5:611138afc79e8da4f1ebda37f6c51ad2</id>
    <published>2420-05-07T15:06:00+01:00</published>
    <updated>2012-05-26T11:51:24+02:00</updated>
    <author><name>thierry sousbie</name></author>
        <dc:subject>Skeleton data</dc:subject>
            
    <content type="html">    &lt;p&gt;VTK formats are developed for the Visualization Tool Kit library (&lt;a href=&quot;http://localhost/dotclear/index.php?post/vtk-formats&quot;&gt;VTK&lt;/a&gt;) and can be used for 3D visualization with software such as &lt;a href=&quot;https://wci.llnl.gov/codes/visit/&quot;&gt;VisIt&lt;/a&gt; or &lt;a href=&quot;http://www.paraview.org/&quot;&gt;ParaView&lt;/a&gt;. Skeletons are stored as &lt;em&gt;VTK polygon data&lt;/em&gt; and can be output in fours different VTK formats:
&lt;br /&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;vtk&lt;/strong&gt;: the legacy format&lt;/li&gt;
&lt;li&gt;-&lt;strong&gt;vtk_ascii&lt;/strong&gt;: ASCII version of the vtk format&lt;/li&gt;
&lt;li&gt;-&lt;strong&gt;vtp&lt;/strong&gt;: a more recently developed XML version of the vtk format,&lt;/li&gt;
&lt;li&gt;-&lt;strong&gt;vtp_ascii&lt;/strong&gt;: ASCII version of the vtp format&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;
&lt;br /&gt;
The specifications for these formats can be found in this &lt;a href=&quot;http://localhost/dotclear/index.php?post/www.vtk.org/VTK/img/file-formats.pdf&quot; title=&quot;VTK formats PDF&quot;&gt;PDF&lt;/a&gt; file. See also &lt;a href=&quot;http://www.cacr.caltech.edu/~slombey/asci/vtk/vtk_formats.simple.html&quot;&gt;here&lt;/a&gt; and &lt;a href=&quot;http://mathema.tician.de/node/430&quot; title=&quot;there&quot;&gt;there&lt;/a&gt; for additional information.
&lt;br /&gt;
&lt;br /&gt;
&lt;br /&gt;
&lt;br /&gt;
&lt;br /&gt;&lt;/p&gt;</content>
    
    

    
      </entry>
    
  <entry>
    <title>crits_ascii format</title>
    <link href="http://localhost/dotclear/index.php?post/crits_ascii-format" rel="alternate" type="text/html"
    title="crits_ascii format" />
    <id>urn:md5:a36d6123b06f743a8cfd629740502e4d</id>
    <published>2420-05-06T15:03:00+01:00</published>
    <updated>2012-05-24T16:31:02+02:00</updated>
    <author><name>thierry sousbie</name></author>
        <dc:subject>Skeleton data</dc:subject>
            
    <content type="html">    &lt;p&gt;This ASCII format is a very simple one. It consists in a list of all the critical points with some useful information for each of them.
&lt;br /&gt;
&lt;br /&gt;&lt;/p&gt;
&lt;table border=&quot;1&quot;&gt;&lt;caption&gt; crits_ascii format&lt;/caption&gt;&lt;tr &gt;&lt;td  width=&quot;50%&quot; &gt;&lt;/td&gt;&lt;td  width=&quot;50%&quot; &gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; #critical points&lt;/th&gt;&lt;td&gt; Header&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; #X0 X1 X2 value type pair_id boundary&lt;/th&gt;&lt;td&gt; Header identifying each column&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; #NDIMS NPOINTS&lt;/th&gt;&lt;td&gt; Number of dimensions and number of points&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; X0 X1 X2 value type pair_id boundary&lt;/th&gt;&lt;td&gt; The value of the fields for the first point&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; ....&lt;/th&gt;&lt;td&gt; One such line for each NPOINTS points&lt;/td&gt;&lt;/tr&gt;
&lt;/table&gt;

&lt;p&gt;&lt;br /&gt;
&lt;br /&gt;
&lt;ins&gt;Notations&lt;/ins&gt;:&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;X_0&lt;/strong&gt;,..., &lt;strong&gt;X_ndims&lt;/strong&gt; : the coordinates of the pooint&lt;/li&gt;
&lt;li&gt;-&lt;strong&gt;value&lt;/strong&gt; : Value of the field.&lt;/li&gt;
&lt;li&gt;-&lt;strong&gt;type&lt;/strong&gt; : the critical index of the point. A value of NDIMS+1 indicates a bifurcation point or a trimmed extremity.&lt;/li&gt;
&lt;li&gt;-&lt;strong&gt;pair_id&lt;/strong&gt; : index of the other critical point in the persistence pair. C style convention is used, indices starting at 0, and points without perisstence pairs have their own index for &lt;em&gt;pair_id&lt;/em&gt;.&lt;/li&gt;
&lt;li&gt;-&lt;strong&gt;boundary&lt;/strong&gt; : value can be &lt;em&gt;0&lt;/em&gt;, &lt;em&gt;1&lt;/em&gt; or &lt;em&gt;2&lt;/em&gt; when the point is inside the domain, on the boundary, or outside (e.g. at infinity)&lt;/li&gt;
&lt;/ul&gt;</content>
    
    

    
      </entry>
    
  <entry>
    <title>segs_ascii format</title>
    <link href="http://localhost/dotclear/index.php?post/segs_ascii-format" rel="alternate" type="text/html"
    title="segs_ascii format" />
    <id>urn:md5:364e501555e6d066eda82436e0a95325</id>
    <published>2420-05-04T15:03:00+01:00</published>
    <updated>2013-01-24T04:56:07+01:00</updated>
    <author><name>thierry sousbie</name></author>
        <dc:subject>Skeleton data</dc:subject>
            
    <content type="html">    &lt;p&gt;This ASCII format is a very simple one. It consists in a list of individual segments describing the local orientation of the arcs as well as some limited information on the filament they belong. This may be used for doing statistical analysis of local properties of filaments such as the local orientation of a magnetic field for instance. Note that all the segments are oriented from the lower critical index extremity of the filament they belong toward the other extremity (which does NOT mean that the value is strictly increasing, as local fluctuations of amplitude the persistence threshold are still allowed ... ).
&lt;br /&gt;
&lt;ins&gt;&lt;strong&gt;nb&lt;/strong&gt;&lt;/ins&gt;: you probably want to use option &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#breakdown&quot;&gt;-breakdown&lt;/a&gt;&lt;/em&gt; of &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv&quot;&gt;skelconv&lt;/a&gt;&lt;/em&gt; before using this format.
&lt;br /&gt;
&lt;br /&gt;&lt;/p&gt;
&lt;table border=&quot;1&quot;&gt;&lt;caption&gt; segs_ascii format&lt;/caption&gt;&lt;tr &gt;&lt;td  width=&quot;50%&quot; &gt;&lt;/td&gt;&lt;td  width=&quot;50%&quot; &gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; #arc segments&lt;/th&gt;&lt;td&gt; Header&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; #U0 U1 U2 V0 V1 V2 value_U value_V type boundary&lt;/th&gt;&lt;td&gt; Header identifying each column&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; #NDIMS NSEGS&lt;/th&gt;&lt;td&gt; Number of dimensions and number of segments&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; U0 U1 U2 V0 V1 V2 value_U value_V type boundary&lt;/th&gt;&lt;td&gt; The value of the fields for the first segment&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; ....&lt;/th&gt;&lt;td&gt; One such line for each NSEGS segment&lt;/td&gt;&lt;/tr&gt;
&lt;/table&gt;

&lt;p&gt;&lt;br /&gt;
&lt;br /&gt;
&lt;ins&gt;Notations&lt;/ins&gt;:&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;U_0&lt;/strong&gt;,..., &lt;strong&gt;U_ndims&lt;/strong&gt; : The coordinates of the first extremity of the segment. It is the one closest to the lowest critical index CP along the filament it belongs to.&lt;/li&gt;
&lt;li&gt;-&lt;strong&gt;V_0&lt;/strong&gt;,..., &lt;strong&gt;V_ndims&lt;/strong&gt; : The coordinates of the second extremity of the segment. It is the one closest to the highest critical index CP along the filament it belongs to.&lt;/li&gt;
&lt;li&gt;-&lt;strong&gt;value_U&lt;/strong&gt;, &lt;strong&gt;value_V&lt;/strong&gt; : Value of the field in U and V.&lt;/li&gt;
&lt;li&gt;-&lt;strong&gt;type&lt;/strong&gt; : the type of arc the segment belongs to. The type of an arc is the minimum critical index of the CP at the extremity of the filament the segment belong to.&lt;/li&gt;
&lt;li&gt;-&lt;strong&gt;boundary&lt;/strong&gt; : value can be &lt;em&gt;0&lt;/em&gt;, &lt;em&gt;1&lt;/em&gt; or &lt;em&gt;2&lt;/em&gt; when the segment is inside the domain, on the boundary, or outside (e.g. at infinity)&lt;/li&gt;
&lt;/ul&gt;</content>
    
    

    
      </entry>
    
  <entry>
    <title>NDskl_ascii format</title>
    <link href="http://localhost/dotclear/index.php?post/NDskl_ascii-format" rel="alternate" type="text/html"
    title="NDskl_ascii format" />
    <id>urn:md5:6f0d89c6789bb8662825d1d4061e3b81</id>
    <published>2420-05-03T15:02:00+01:00</published>
    <updated>2012-05-27T15:50:00+02:00</updated>
    <author><name>thierry sousbie</name></author>
        <dc:subject>Skeleton data</dc:subject>
            
    <content type="html">    &lt;p&gt;An ASCII format that contains the same amount of information as &lt;a href=&quot;http://localhost/dotclear/index.php?post/NDskl-format&quot;&gt;NDskl&lt;/a&gt; files, but organized in a different way. In particular, filaments are not described as lists of segments but rather each filament is described by an origin node, a destination node, and a set of sampling points.
&lt;br /&gt;
&lt;br /&gt;&lt;/p&gt;
&lt;table border=&quot;1&quot;&gt;&lt;caption&gt; NDskl_ascii format&lt;/caption&gt;&lt;tr &gt;&lt;td  width=&quot;50%&quot; &gt;&lt;/td&gt;&lt;td  width=&quot;50%&quot; &gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; ANDSKEL&lt;/th&gt;&lt;td&gt; header&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; ndims&lt;/th&gt;&lt;td&gt; the number of dimensions&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; #comments go here&lt;/th&gt;&lt;td&gt; OPTIONAL: should start with '#' if present (the 80 first characters are read and stored).&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; BBOX [x0_1 .. x0_d] [delta_1 .. delta_d]&lt;/th&gt;&lt;td&gt; OPTIONAL: the bounding box, defined by the 'ndims' coordinates of the origin 'x0' and extent 'delta'.&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; [CRITICAL POINTS]&lt;/th&gt;&lt;td&gt; Marks the beginning of the critical points section&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; ncrit&lt;/th&gt;&lt;td&gt; The number of critical points (CP)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; type pos[0] ... pos[ndims-1] value pairID boundary&lt;/th&gt;&lt;td style=&quot;background:lightblue&quot;&gt; Info on the first CP: critical index, position, value, index of CP in the persistence pair, 0 if not on the boundary&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; nfil&lt;/th&gt;&lt;td style=&quot;background:lightblue&quot;&gt; The number of filaments connected to this CP&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; destId[0] filId[0]&lt;/th&gt;&lt;td style=&quot;background:lightblue&quot;&gt; Info on the first filament: index of the CP at the other extremity of the filament, and index of the filament (see filaments table below)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; ...&lt;/th&gt;&lt;td style=&quot;background:lightblue&quot;&gt; One line for each filament connecting on the CP&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; destId[nfil-1] filId[nfil-1]&lt;/th&gt;&lt;td style=&quot;background:lightblue&quot;&gt; Information on the last filament&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th  style=&quot;background:red&quot;&gt;.....&lt;/th&gt;&lt;td  style=&quot;background:red&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th style=&quot;background:red&quot;&gt; .....&lt;/th&gt;&lt;td style=&quot;background:red&quot;&gt; one blue block for each CP.&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th style=&quot;background:red&quot;&gt; .....&lt;/th&gt;&lt;td style=&quot;background:red&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; [FILAMENTS]&lt;/th&gt;&lt;td&gt; Marks the beginning of the filaments section&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; nfil&lt;/th&gt;&lt;td&gt; Total number of filaments&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; CP1 CP2 nSamp&lt;/th&gt;&lt;td style=&quot;background:lightblue&quot;&gt; index of the CP at the extremity of the first filament and number of sampling points&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; P[0][0] ... P[0][ndims-1]&lt;/th&gt;&lt;td  style=&quot;background:lightblue&quot;&gt;position of the first sampling point of first filament.&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; ...&lt;/th&gt;&lt;td  style=&quot;background:lightblue&quot;&gt;One line for each sampling point of first filament.&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; P[nSamp-1][0] ... P[nSamp-1][ndims-1]&lt;/th&gt;&lt;td  style=&quot;background:lightblue&quot;&gt;Position of the last sampling point&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th  style=&quot;background:red&quot;&gt;.....&lt;/th&gt;&lt;td  style=&quot;background:red&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th style=&quot;background:red&quot;&gt; .....&lt;/th&gt;&lt;td style=&quot;background:red&quot;&gt; one blue block for each filament.&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th style=&quot;background:red&quot;&gt; .....&lt;/th&gt;&lt;td style=&quot;background:red&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; [CRITICAL POINTS DATA]&lt;/th&gt;&lt;td&gt; Marks the beginning of the CP data section&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; NF&lt;/th&gt;&lt;td&gt; Number of fields associated to each CP.&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; CP_DATA_FIELD_NAME_1&lt;/th&gt;&lt;td&gt; Name of the first field&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; ...&lt;/th&gt;&lt;td&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; CP_DATA_FIELD_NAME_NF&lt;/th&gt;&lt;td&gt; Name of the last field&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; val_1[0] ... val_NF[0]&lt;/th&gt;&lt;td&gt; Value of each field for first CP&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; ...&lt;/th&gt;&lt;td&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; val_1[N_CP-1] ... val_NF[N_CP-1]&lt;/th&gt;&lt;td&gt; Value of each field for last CP&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; [FILAMENTS DATA]&lt;/th&gt;&lt;td&gt; Marks the beginning of the filaments data section&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; NF&lt;/th&gt;&lt;td&gt; Number of fields associated to each sampling point of each filament.&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; FIL_DATA_FIELD_NAME_1&lt;/th&gt;&lt;td&gt; Name of the first field&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; ...&lt;/th&gt;&lt;td&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; FIL_DATA_FIELD_NAME_NF&lt;/th&gt;&lt;td&gt; Name of the last field&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; val_1[0][0] ... val_NF[0][0]&lt;/th&gt;&lt;td style=&quot;background:lightblue&quot;&gt; field values for first sampling point, first filament&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; ...&lt;/th&gt;&lt;td style=&quot;background:lightblue&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; val_1[0][nSamp[0]] ... val_NF[0][nSamp[0]]&lt;/th&gt;&lt;td style=&quot;background:lightblue&quot;&gt; field values for last sampling point, first filament&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th  style=&quot;background:red&quot;&gt;.....&lt;/th&gt;&lt;td  style=&quot;background:red&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th style=&quot;background:red&quot;&gt; .....&lt;/th&gt;&lt;td style=&quot;background:red&quot;&gt; one blue block for each filament.&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th style=&quot;background:red&quot;&gt; .....&lt;/th&gt;&lt;td style=&quot;background:red&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;/table&gt;</content>
    
    

    
      </entry>
    
  <entry>
    <title>NDskl format</title>
    <link href="http://localhost/dotclear/index.php?post/NDskl-format" rel="alternate" type="text/html"
    title="NDskl format" />
    <id>urn:md5:c09ea460592168c16755715666249515</id>
    <published>2420-05-02T15:01:00+01:00</published>
    <updated>2012-05-27T12:55:36+02:00</updated>
    <author><name>thierry sousbie</name></author>
        <dc:subject>Skeleton data</dc:subject>
            
    <content type="html">    &lt;p&gt;This is the native binary format of DisPerSE. Functions for reading and writing &lt;em&gt;NDskl&lt;/em&gt; format  in &lt;em&gt;C&lt;/em&gt; can be found within the file &lt;code&gt;${DISPERSE_SRC}/src/C/NDskeleton.c&lt;/code&gt; (see functions &lt;em&gt;Load_NDskel&lt;/em&gt;, &lt;em&gt;Save_NDskel&lt;/em&gt; and also &lt;em&gt;Create_NDskel&lt;/em&gt;). In this format, the arcs of the Morse-Smale complex are described as arrays of nodes and segments with data associated to each of them. Each node contains information on critical points (which may also be bifurcations points or trimmed extremities of filaments). Each segment contains information on the local geometry of the arcs and global properties of the skeleton.
&lt;br /&gt;
&lt;br /&gt;
When using the C functions from Disperse, data is loaded into the following &lt;em&gt;C&lt;/em&gt; structure which is close to the actual structure of the file (see file &lt;code&gt;${DISPERSE_SRC}/src/C/NDskeleton.h&lt;/code&gt;):&lt;/p&gt;


&lt;pre&gt;struct NDskl_seg_str {
   int nodes[2];  // index of the nodes at the extremity of the arc. Segment is oriented from nodes[0] toward nodes[1]
   float *pos;  // points to appropriate location in segpos
   int flags;  // non null if on boundary 
   int index;  // index of the segment in the Seg array
   double *data;  // points to the nsegdata supplementary data for this segment
   struct NDskl_seg_str *Next; // points to next segment in the arc, NULL if extremity (-&amp;gt;connects to nodes[1])
   struct NDskl_seg_str *Prev; // points to previous segment in the arc, NULL if extremity (-&amp;gt;connects to nodes[0])
 }; 
 typedef struct NDskl_seg_str NDskl_seg;&lt;/pre&gt;


&lt;pre&gt; struct NDskl_node_str {
   float *pos;  // points to appropriate location in nodepos
   int flags;  // non null if on boundary 
   int nnext;  // number of arcs connected
   int type;  // critical index
   int index;  // index of the node in the Node array
   int *nsegs; // number pf segments in each of the the nnext arcs
   double *data;   // points to the nnodedata supplementary data for this segment
   struct NDskl_node_str **Next;  // points to the node at the other extremity of each of the nnext arcs
   struct NDskl_seg_str **Seg;   // points to the first segment in the arcs, starting from current node.
 }; 
 typedef struct NDskl_node_str NDskl_node;&lt;/pre&gt;


&lt;pre&gt; typedef struct NDskel_str {
   char comment[80];
   
   int ndims;  // number of dimensions
   int *dims;  // dimension of the underlying grid, only meaningfull when computed from a regular grid
   double *x0;  // origin of the bounding box
   double *delta;  // dimension of the bounding box
    
   int nsegs;  // number of segments
   int nnodes;  // number of nodes
   
   int nsegdata;  // number of additional data for segments
   int nnodedata;  // number of additional data for nodes
   char **segdata_info;  // name of additional fields for segments
   char **nodedata_info;  // name of additional fields for nodes
   
   float *segpos;  // positions of the extremities of segments (2xndims coords for each segment)
   float *nodepos; // positions of the nodes (ndims coords for each segment)
   double *segdata; // additional data for segments (nsegs times nsegdata consecutive values)
   double *nodedata; // additional data for nodes (nnodes times nnodedata consecutive values)
   
   NDskl_seg *Seg;  // segment array (contains all segs)
   NDskl_node *Node;  // nodes array (contains all nodes)
 } NDskel;&lt;/pre&gt;

&lt;p&gt;&lt;br /&gt;
&lt;br /&gt;
The &lt;em&gt;NDskl&lt;/em&gt; binary format is organized as follows (blocks are delimited by &lt;em&gt;dummy&lt;/em&gt; variables indicating the size of the blocks for FORTRAN compatibility, but they are ignored in C):
&lt;br /&gt;
&lt;br /&gt;&lt;/p&gt;
&lt;table style=&quot;margin: 1em auto 1em auto&quot;&gt;&lt;caption&gt;NDskl binary format&lt;/caption&gt;&lt;tr  valign=&quot;top&quot;&gt;&lt;td width=&quot;10%&quot;&gt;field&lt;/td&gt;&lt;td width=&quot;10%&quot;&gt;type&lt;/td&gt;&lt;td width=&quot;10%&quot;&gt;size&lt;/td&gt;&lt;td width=&quot;65%&quot;&gt;comment&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt; for FORTRAN compatibility&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;tag&lt;/td&gt;&lt;td&gt;char(1B)&lt;/td&gt;&lt;td&gt;16&lt;/td&gt;&lt;td&gt;identifies the file type. Value : &quot;NDSKEL&quot;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;comment&lt;/td&gt;&lt;td&gt;char(1B)&lt;/td&gt;&lt;td&gt;80&lt;/td&gt;&lt;td&gt;a comment on the file (string)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;ndims&lt;/td&gt;&lt;td&gt;int(4B)&lt;/td&gt;&lt;td&gt;1&lt;/td&gt;&lt;td&gt;number of dimensions&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;dims&lt;/td&gt;&lt;td&gt;int(4B)&lt;/td&gt;&lt;td&gt;20&lt;/td&gt;&lt;td&gt;dimension of underlying grid (0=none, ndims values are read)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;x0&lt;/td&gt;&lt;td&gt;double(8B)&lt;/td&gt;&lt;td&gt;20&lt;/td&gt;&lt;td&gt;origin of bounding box (ndims first values are meaningful)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;delta&lt;/td&gt;&lt;td&gt;double(8B)&lt;/td&gt;&lt;td&gt;20&lt;/td&gt;&lt;td&gt;size of bounding box (ndims first values are meaningful)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;nsegs&lt;/td&gt;&lt;td&gt;int(4B)&lt;/td&gt;&lt;td&gt;1&lt;/td&gt;&lt;td&gt;number of segments&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;nnodes&lt;/td&gt;&lt;td&gt;int(4B)&lt;/td&gt;&lt;td&gt;1&lt;/td&gt;&lt;td&gt;number of nodes&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;nsegdata&lt;/td&gt;&lt;td&gt;int(4B)&lt;/td&gt;&lt;td&gt;1&lt;/td&gt;&lt;td&gt;number of additional data associated to segments.&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;nnodedata&lt;/td&gt;&lt;td&gt;int(4B)&lt;/td&gt;&lt;td&gt;1&lt;/td&gt;&lt;td&gt;number of additional data associated to nodes.&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;  this block is ommited if nsegdata=0&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;seg_data_info&lt;/td&gt;&lt;td&gt;char(1B)&lt;/td&gt;&lt;td&gt;20&amp;timesnsegdata&lt;/td&gt;&lt;td&gt;name of the segment data field&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt; this block is ommited if nsegdata=0&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;  this block is ommited if nnodedata=0&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;node_data_info&lt;/td&gt;&lt;td&gt;char(1B)&lt;/td&gt;&lt;td&gt;20&amp;timesnnodedata&lt;/td&gt;&lt;td&gt;name of the node data field&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;  this block is ommited if nnodedata=0&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;segpos&lt;/td&gt;&lt;td&gt;float(4B)&lt;/td&gt;&lt;td&gt;2&amp;timesnsegs&amp;timesndims&lt;/td&gt;&lt;td&gt;coordinates of segment extremities&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;nodepos&lt;/td&gt;&lt;td&gt;float(4B)&lt;/td&gt;&lt;td&gt;nnodes&amp;timesndims&lt;/td&gt;&lt;td&gt;coordinates of nodes&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;segdata&lt;/td&gt;&lt;td&gt;double(8B)&lt;/td&gt;&lt;td&gt;nsegdata&amp;timesnsegs&lt;/td&gt;&lt;td&gt;value of the segments data fields&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;nodedata&lt;/td&gt;&lt;td&gt;double(8B)&lt;/td&gt;&lt;td&gt;nnodedata&amp;timesnnodes&lt;/td&gt;&lt;td&gt;value of the nodes data fields&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:red&quot;&gt; following block is repeated for each node (&amp;timesnnodes).&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;pos_index&lt;/td&gt;&lt;td&gt;int(4B)&lt;/td&gt;&lt;td&gt;1&lt;/td&gt;&lt;td&gt;index in the nodepos/nodedata arrays&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;flags&lt;/td&gt;&lt;td&gt;int(4B)&lt;/td&gt;&lt;td&gt;1&lt;/td&gt;&lt;td&gt;flags (identify boundary nodes)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;nnext&lt;/td&gt;&lt;td&gt;int(4B)&lt;/td&gt;&lt;td&gt;1&lt;/td&gt;&lt;td&gt;number of connected arcs&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;type&lt;/td&gt;&lt;td&gt;int(4B)&lt;/td&gt;&lt;td&gt;1&lt;/td&gt;&lt;td&gt;critical index (ndims+1 for bifurcations / trimmed arc axtremity)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;index&lt;/td&gt;&lt;td&gt;int(4B)&lt;/td&gt;&lt;td&gt;1&lt;/td&gt;&lt;td&gt;index of this node&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;nsegs&lt;/td&gt;&lt;td&gt;int(4B)&lt;/td&gt;&lt;td&gt;nnext&lt;/td&gt;&lt;td&gt;number of segments in each connected arc&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt; * &lt;/td&gt;&lt;td&gt;*&lt;/td&gt;&lt;td&gt;*&lt;/td&gt;&lt;td  style=&quot;background:red&quot;&gt;blue section repeated for each connected arc (&amp;timesnnext)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;nextNode&lt;/td&gt;&lt;td&gt;int(4B)&lt;/td&gt;&lt;td&gt;1&lt;/td&gt;&lt;td style=&quot;background:lightblue&quot;&gt;index of the other node of the arc&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;nextSeg&lt;/td&gt;&lt;td&gt;int(4B)&lt;/td&gt;&lt;td&gt;1&lt;/td&gt;&lt;td style=&quot;background:lightblue&quot;&gt;index of the first segment on the arc&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:red&quot;&gt; following block is repeated for each segment (&amp;timesnsegs).&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;pos_index&lt;/td&gt;&lt;td&gt;int(4B)&lt;/td&gt;&lt;td&gt;1&lt;/td&gt;&lt;td&gt;index in the segpos/segdata arrays&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;nodes&lt;/td&gt;&lt;td&gt;int(4B)&lt;/td&gt;&lt;td&gt;2&lt;/td&gt;&lt;td&gt;index of the 2 nodes at the endpoints of the current arc.&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;flags&lt;/td&gt;&lt;td&gt;int(4B)&lt;/td&gt;&lt;td&gt;1&lt;/td&gt;&lt;td&gt;flags (identify boundary nodes)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;index&lt;/td&gt;&lt;td&gt;int(4B)&lt;/td&gt;&lt;td&gt;1&lt;/td&gt;&lt;td&gt;index of this segment&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;next_seg&lt;/td&gt;&lt;td&gt;int(4B)&lt;/td&gt;&lt;td&gt;1&lt;/td&gt;&lt;td&gt;index of the next segment in the node (-1 if none)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;prev_seg&lt;/td&gt;&lt;td&gt;int(4B)&lt;/td&gt;&lt;td&gt;1&lt;/td&gt;&lt;td&gt;index of the previous segment in the node (-1 if none)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;/tr&gt;
&lt;/table&gt;</content>
    
    

    
      </entry>
    
  <entry>
    <title>Skeleton files</title>
    <link href="http://localhost/dotclear/index.php?post/skeleton-formats" rel="alternate" type="text/html"
    title="Skeleton files" />
    <id>urn:md5:5a11e7d49bef99b40d28c01c537d6cad</id>
    <published>2420-05-01T06:56:00+01:00</published>
    <updated>2012-07-18T11:10:05+02:00</updated>
    <author><name>thierry sousbie</name></author>
        <dc:subject>Skeleton data</dc:subject>
            
    <content type="html">    &lt;p&gt;This file type is designed to store the critical points and arcs of the Morse-Smale complex (i.e. one dimensionnal filamentary structures). In skeleton files, the filamentary network is described as lists of nodes representing critical points or bifurcation points if available (see option &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#breakdown&quot;&gt;-breakdown&lt;/a&gt;&lt;/em&gt; of &lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv&quot;&gt;skelconv&lt;/a&gt;), lists of segments tracing the local geometry of arcs and filaments originating from and leading to nodes and described as lists of segments. The base skeleton format is &lt;a href=&quot;http://localhost/dotclear/index.php?post/NDskl-format&quot;&gt;NDskl&lt;/a&gt; which is used internally, but this format may be converted to several other more or less complex formats of skeleton files adapted to different applications (see option &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#to&quot;&gt;-to&lt;/a&gt;&lt;/em&gt; in program &lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv&quot;&gt;skelconv&lt;/a&gt;, a list of available formats is displayed when running the program without argument).&lt;br /&gt;
&lt;ins&gt;nb&lt;/ins&gt;: Within skeleton files, the segments are always oriented in the ascending direction of the arcs, which does *NOT* necessarily mean that the value of the field is necessarily increasing locally ! Indeed, the value is locally allowed to fluctuate within the persistence threshold, so the field value does not have to be strictly increasing locally along an arc.
&lt;br /&gt;
&lt;br /&gt;
&lt;strong&gt;&lt;ins&gt;Available formats&lt;/ins&gt;&lt;/strong&gt;:
&lt;br /&gt;
&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/NDskl-format&quot;&gt;NDskl&lt;/a&gt;&lt;/strong&gt; (Read / Write):&lt;br /&gt; This is the format of the skeleton files created with &lt;a href=&quot;http://localhost/dotclear/index.php?post/mse&quot;&gt;mse&lt;/a&gt;. It is a relatively complex binary format that contains all the information on the geometry and connectivity of the arcs of the Morse-Smale complex (i.e. filaments).&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/NDskl_ascii-format&quot;&gt;NDskl_ascii&lt;/a&gt;&lt;/strong&gt; (Write only):&lt;br /&gt; This ASCII format contains the same amount of information as &lt;em&gt;NDskl&lt;/em&gt; files, but organized in a different way. In particular, filaments are not described as lists of segments but rather each filament is described by an origin node, a destination node, and a set of sampling points.&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/segs_ascii-format&quot;&gt;segs_ascii&lt;/a&gt;&lt;/strong&gt; and &lt;strong&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/crits_ascii-format&quot;&gt;crits_ascii&lt;/a&gt;&lt;/strong&gt; (Write only):&lt;br /&gt; This is a simplified ASCII file format designed to be easily readable but that only contains a subset of the information available in &lt;em&gt;NDskl&lt;/em&gt; files. In particular, &lt;em&gt;segs_ascii&lt;/em&gt; files contain a list of individual segments describing the local orientation of the arcs as well as some limited information on the filament they belong to, while &lt;em&gt;crits_ascii&lt;/em&gt; contain information on the critical points and bifurcation points only.&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/vtk-formats&quot;&gt;vtk&lt;/a&gt;&lt;/strong&gt;, &lt;strong&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/vtk-formats&quot;&gt;vtk_ascii&lt;/a&gt;&lt;/strong&gt;, &lt;strong&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/vtk-formats&quot;&gt;vtp&lt;/a&gt;&lt;/strong&gt; and &lt;strong&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/vtk-formats&quot;&gt;vtp_ascii&lt;/a&gt;&lt;/strong&gt; (Write only):&lt;br /&gt; These formats are binary and ASCII legacy and XML &lt;a href=&quot;http://www.vtk.org/&quot;&gt;VTK&lt;/a&gt; formats that are readable by several 3D visualization tools, such as &lt;a href=&quot;https://wci.llnl.gov/codes/visit/&quot;&gt;VisIt&lt;/a&gt; or &lt;a href=&quot;http://www.paraview.org/&quot;&gt;ParaView&lt;/a&gt; for instance.&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/networks&quot;&gt;NDnet&lt;/a&gt;&lt;/strong&gt; (Write only):&lt;br /&gt; Using this format, skeleton files can be converted to &lt;a href=&quot;http://localhost/dotclear/index.php?post/networks&quot;&gt;unstructured networks&lt;/a&gt; (use &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/netconv&quot;&gt;netconv&lt;/a&gt;&lt;/em&gt; to convert &lt;em&gt;NDnet&lt;/em&gt; files to other network formats).&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;
&lt;br /&gt;
&lt;strong&gt;&lt;ins&gt;Additional data&lt;/ins&gt;&lt;/strong&gt;: In addition to the topology and geometry of filamentary structures, more complete formats (i.e. &lt;a href=&quot;http://localhost/dotclear/index.php?post/NDskl-format&quot;&gt;NDskl&lt;/a&gt;, &lt;a href=&quot;http://localhost/dotclear/index.php?post/NDskl_ascii-format&quot;&gt;NDskl_ascii&lt;/a&gt; and all &lt;a href=&quot;http://localhost/dotclear/index.php?post/vtk-formats&quot;&gt;vtk&lt;/a&gt; formats) may store arbitrary additional information associated to filaments and nodes that can be used by &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv&quot;&gt;skelconv&lt;/a&gt;&lt;/em&gt; (run &lt;em&gt;skelconv filename.NDskl -info&lt;/em&gt; for a list of additional data available in skeleton file &lt;em&gt;filename.NDskl&lt;/em&gt;). By default, &lt;a href=&quot;http://localhost/dotclear/index.php?post/mse&quot;&gt;mse&lt;/a&gt; stores the following additional data in skeleton type files:
&lt;br /&gt;
&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;persistence&lt;/strong&gt; / &lt;strong&gt;persistence_ratio&lt;/strong&gt; / &lt;strong&gt;persistence_nsigmas&lt;/strong&gt; (Nodes only) :&lt;br /&gt; The persistence (expressed as a difference, ratio or in &lt;em&gt;number of sigmas&lt;/em&gt;) of the persistence pair containing the corresponding critical point. A negative or null value indicates that the node is not a critical point or that it does not belong to a persistence pair.&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;persistence_pair&lt;/strong&gt; (Nodes only):&lt;br /&gt; The index of the node that corresponds to the other critical point in the persistence pair (indices start at 0). The value is the index of the current node itself when it is not a critical point or it does not belong to a persistence pair.&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;parent&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;parent_index&lt;/strong&gt; / &lt;strong&gt;parent_log_index&lt;/strong&gt; (Nodes only):&lt;br /&gt; For each extrema only (i.e. minima and maxima), the index of the node that corresponds to the other extremum into which it would be merged if its persistence pair was canceled (indices start at 0). This can be used to reconstruct the tree of the hierarchy of maxima and minima. The value is -1 for non extrema critical points. The difference between the two versions is that the second (&lt;em&gt;parent_log_index&lt;/em&gt;) is the hierarchy computed from the logarithm of the field. This second version is useful only for discrete point samples whose MS-complex is obtained from the delaunay tessellation computed with &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/Usage&quot;&gt;delaunay_nD&lt;/a&gt;&lt;/em&gt;. Practically, &lt;em&gt;parent_log_index&lt;/em&gt; can be used whenever persistence pairs are cancelled in order of increasing ratio (option &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#nsig&quot;&gt;-nsig&lt;/a&gt;&lt;/em&gt; in &lt;em&gt;mse&lt;/em&gt;), and &lt;em&gt;parent_index&lt;/em&gt; whenever persistence pairs are cancelled in order of increasing difference  (option &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#cut&quot;&gt;-cut&lt;/a&gt;&lt;/em&gt; in &lt;em&gt;mse&lt;/em&gt;).&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;field_value&lt;/strong&gt; / &lt;strong&gt;log_field_value&lt;/strong&gt; (Nodes and Segments) :&lt;br /&gt; The value of the field and it logarithm. For segments, the suffix &lt;em&gt;_p1&lt;/em&gt; and &lt;em&gt;_p2&lt;/em&gt; is added to indicate which extremity of the segment it corresponds to.&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;cell&lt;/strong&gt; (Nodes and Segments) :&lt;br /&gt; The &lt;em&gt;type&lt;/em&gt; and &lt;em&gt;index&lt;/em&gt; of the cell corresponding to the node / segment in the initial network (i.e. from which the skeleton was computed). For segments, the suffix &lt;em&gt;_p1&lt;/em&gt; and &lt;em&gt;_p2&lt;/em&gt; is added to indicate which extremity of the segment it corresponds to. The value is a double precision floating number whose integer part is the index of the cell and decimal part its type. For instance, the 156th vertex (i.e. 0-cell) in the cell complex is represented as 156.0, while the 123th tetrahedron is 123.3. Note that the index of the 0-cell correpond to the index of the pixel / vertices in the original network from which the skeleton was computed.&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;robustness&lt;/strong&gt; / &lt;strong&gt;robustness_ratio&lt;/strong&gt; (Nodes and Segments) :&lt;br /&gt; The robustness of the node / segment. Robustness is a local measure of how contrasted the critical point / filament is with respect to its &lt;em&gt;local&lt;/em&gt; background, and it is therefore a good way to select only highly contrasted subsets of the filamentary structures (see option &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv#trim&quot;&gt;-trimBelow&lt;/a&gt;&lt;/em&gt;  in &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv&quot;&gt;skelconv&lt;/a&gt;&lt;/em&gt;). Note that robustness, like persistence, is defined as a difference or ratio between the value of the field in two points, so the robustness threshold has the same order of magnitude as the persistence threshold.&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;- &lt;strong&gt;type&lt;/strong&gt; (Segments only) :&lt;br /&gt; The type of arc the segment belongs to. The value corresponds to the lowest critical index of the two critical points at the extremities of the arc the segments belongs to. For instance, in dimension D, the ridges (filaments) have type &lt;em&gt;D-1&lt;/em&gt;.&lt;/li&gt;
&lt;/ul&gt;</content>
    
    

    
      </entry>
  
</feed>