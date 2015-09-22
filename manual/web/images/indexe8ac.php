<?xml version="1.0" encoding="utf-8"?><feed xmlns="http://www.w3.org/2005/Atom"
  xmlns:dc="http://purl.org/dc/elements/1.1/"
  xmlns:wfw="http://wellformedweb.org/CommentAPI/"
  xml:lang="en">
  
  <title type="html">DisPerSE - persistent structures identification - Network data</title>
  <subtitle type="html">Automatic identification of persistent structures in 2D or 3D.
keywords: Morse complex, topology, peak, void, source, wall, filament, cosmic web, cosmology, structure identification.</subtitle>
  <link href="http://localhost/dotclear/index.php?feed/category/Network-I-O/atom" rel="self" type="application/atom+xml"/>
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
    <title>vtk network format</title>
    <link href="http://localhost/dotclear/index.php?post/vtk-network-format" rel="alternate" type="text/html"
    title="vtk network format" />
    <id>urn:md5:fb5abf70ccbb5abba106e974c9f74000</id>
    <published>2426-01-04T13:14:00+01:00</published>
    <updated>2012-05-27T19:26:02+02:00</updated>
    <author><name>thierry sousbie</name></author>
        <dc:subject>Network data</dc:subject>
            
    <content type="html">    &lt;p&gt;VTK formats are developed for the Visualization Tool Kit library (&lt;a href=&quot;http://localhost/dotclear/index.php?post/vtk-formats&quot;&gt;VTK&lt;/a&gt;) and can be used for 3D visualization with software such as &lt;a href=&quot;https://wci.llnl.gov/codes/visit/&quot;&gt;VisIt&lt;/a&gt; or &lt;a href=&quot;http://www.paraview.org/&quot;&gt;ParaView&lt;/a&gt;. Networks are stored as &lt;em&gt;VTK unstructured network data&lt;/em&gt; and can be output in fours different VTK formats:
&lt;br /&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;vtk&lt;/strong&gt;: the legacy format&lt;/li&gt;
&lt;li&gt;-&lt;strong&gt;vtk_ascii&lt;/strong&gt;: ASCII version of the vtk format&lt;/li&gt;
&lt;li&gt;-&lt;strong&gt;vtu&lt;/strong&gt;: a more recently developed XML version of the vtk format,&lt;/li&gt;
&lt;li&gt;-&lt;strong&gt;vtu_ascii&lt;/strong&gt;: ASCII version of the vtu format&lt;/li&gt;
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
    <title>ply format</title>
    <link href="http://localhost/dotclear/index.php?post/ply-format" rel="alternate" type="text/html"
    title="ply format" />
    <id>urn:md5:776d9154fe52b9084d1ffcd74bdebf79</id>
    <published>2426-01-03T13:14:00+01:00</published>
    <updated>2012-05-27T19:27:27+02:00</updated>
    <author><name>thierry sousbie</name></author>
        <dc:subject>Network data</dc:subject>
            
    <content type="html">    &lt;p&gt;The PLY format is a relatively generic file format designed to store three dimensional data from 3D scanners with the possibility of associating properties to the polygons. Information on this format can be found on &lt;a href=&quot;http://en.wikipedia.org/wiki/PLY_(file_format)&quot; title=&quot;PLY format&quot;&gt;wikipedia&lt;/a&gt;  (see also the &lt;em&gt;External links&lt;/em&gt; section). A very efficient C library for reading and writing PLY files in ASCII or binary format is &lt;a href=&quot;http://w3.impa.br/~diego/software/rply/&quot; title=&quot;RPly&quot;&gt;RPly&lt;/a&gt; (DisPerSE uses it for PLY files I/O, see files &lt;em&gt;NDnet_PLY_IO.c&lt;/em&gt; and &lt;em&gt;NDnet_PLY_IO.h&lt;/em&gt;).
&lt;br /&gt;
&lt;br /&gt;
A typical header for a PLY file readable by &lt;a href=&quot;http://localhost/dotclear/index.php?post/netconv&quot;&gt;netconv&lt;/a&gt; or &lt;a href=&quot;http://localhost/dotclear/index.php?post/mse&quot;&gt;mse&lt;/a&gt; is as follows:&lt;/p&gt;


&lt;pre&gt;ply
format ascii 1.0
element bbox 1
property list uchar double x0 
property list uchar double delta
element vertex 77595
property float x
property float y
property float z
property double field_value
element face 482867
property list uchar uint vertex_indices
end_header&lt;/pre&gt;

&lt;p&gt;&lt;br /&gt;
&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;format&lt;/strong&gt; may be ASCII or little / big endian binary&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;bbox&lt;/strong&gt; element is used to define a bounding box if available (x0 is its origin and delta its extent)&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;coordinates&lt;/strong&gt; of the vertices are given as vertex properties labeled &lt;em&gt;x&lt;/em&gt;, &lt;em&gt;y&lt;/em&gt;, and &lt;em&gt;z&lt;/em&gt; or &lt;em&gt;x0&lt;/em&gt;, &lt;em&gt;x1&lt;/em&gt;, ... (the number of this properties gives the dimension of the embedding space)&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;faces&lt;/strong&gt; are defined by the property vertex_indices, each element corresponding to a list of vertices. A cell is always supposed to be a simplex, so the number of vertices determine the dimension spanned by the complex (a 2D complex may be embedded in a 3D space).&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;additional properties&lt;/strong&gt; may be defined for cells and vertices. In particular, a  vertex property labeled &lt;em&gt;field_value&lt;/em&gt; would be used as input function in &lt;a href=&quot;http://localhost/dotclear/index.php?post/mse&quot;&gt;mse&lt;/a&gt;.&lt;/li&gt;
&lt;/ul&gt;</content>
    
    

    
      </entry>
    
  <entry>
    <title>NDnet_ascii format</title>
    <link href="http://localhost/dotclear/index.php?post/NDnet_ascii-format" rel="alternate" type="text/html"
    title="NDnet_ascii format" />
    <id>urn:md5:e274f99d05092ef3465a29a5b43d7631</id>
    <published>2426-01-02T13:14:00+01:00</published>
    <updated>2012-07-18T11:17:00+02:00</updated>
    <author><name>thierry sousbie</name></author>
        <dc:subject>Network data</dc:subject>
            
    <content type="html">    &lt;p&gt;This ASCII format is a simpler version of the &lt;a href=&quot;http://localhost/dotclear/index.php?post/NDnet-format&quot;&gt;NDnet&lt;/a&gt; format, designed to be fully compatible but restricted to simplicial networks. It is easy to read and write and should probably be used for reasonably sized data sets.
&lt;br /&gt;
&lt;br /&gt;
&lt;ins&gt;Note&lt;/ins&gt;: The scalar function whose MS-complex is computed by &lt;a href=&quot;http://localhost/dotclear/index.php?post/mse&quot;&gt;mse&lt;/a&gt; can be stored as an additional data field
named 'field_value' (case sensitive).&lt;/p&gt;


&lt;p&gt;&lt;br /&gt;
&lt;br /&gt;&lt;/p&gt;
&lt;table border=&quot;1&quot;&gt;&lt;caption&gt; NDnet_ascii format&lt;/caption&gt;&lt;tr &gt;&lt;td  width=&quot;50%&quot; &gt;&lt;/td&gt;&lt;td  width=&quot;50%&quot; &gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; ANDNET&lt;/th&gt;&lt;td&gt; header&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; ndims&lt;/th&gt;&lt;td&gt; the number of dimensions&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; #comments go here&lt;/th&gt;&lt;td&gt; OPTIONAL: should start with '#' if present (the 80 first characters are read and stored).&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; BBOX [x0_1 .. x0_d] [delta_1 .. delta_d]&lt;/th&gt;&lt;td&gt; OPTIONAL: the bounding box, defined by the 'ndims' coordinates of the origin 'x0' and extent 'delta'.&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; nv&lt;/th&gt;&lt;td&gt; number of vertices&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; vx[0] vy[0] ...&lt;/th&gt;&lt;td&gt; the ndims coordinates of the first vertex&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; ...&lt;/th&gt;&lt;td&gt; One line for each vertex&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:red&quot; colspan=&quot;2&quot; &gt;Simplices definition (0-cells up to ndims-cells). One blue block should be added for each type of explicitly defined simplex. Note that only the highest dimension cells are sufficient to define a complex.&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th style=&quot;background:lightblue&quot;&gt; T N&lt;/th&gt;&lt;td style=&quot;background:lightblue&quot;&gt; network has N T-simplices (each T-simplex has (T+1) vertices).&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th style=&quot;background:lightblue&quot;&gt; i[0] j[0] ...&lt;/th&gt;&lt;td style=&quot;background:lightblue&quot;&gt; the T+1 indices (start at 0) of the vertices of the first T-simplex.&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th style=&quot;background:lightblue&quot;&gt; ...&lt;/th&gt;&lt;td style=&quot;background:lightblue&quot;&gt; one line for each of the N T-simplices&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; [ADDITIONAL_DATA]&lt;/th&gt;&lt;td&gt; OPTIONAL: indicate the beginning of the additional data section&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; additional_data_name_1&lt;/th&gt;&lt;td&gt; name of the additional data (e.g. field_value for mse input files)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; T&lt;/th&gt;&lt;td&gt; type of simplex it is associated to (T-simplex, 0 means vertices)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; val[0]&lt;/th&gt;&lt;td&gt; value for the first T-simplex&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt; ...&lt;/th&gt;&lt;td&gt; One line for each T-simplex&lt;/td&gt;&lt;/tr&gt;
&lt;/table&gt;</content>
    
    

    
      </entry>
    
  <entry>
    <title>NDnet format</title>
    <link href="http://localhost/dotclear/index.php?post/NDnet-format" rel="alternate" type="text/html"
    title="NDnet format" />
    <id>urn:md5:5ef6902f2833bcc58788d51806dd11ba</id>
    <published>2426-01-01T13:14:00+01:00</published>
    <updated>2013-01-24T04:57:20+01:00</updated>
    <author><name>thierry sousbie</name></author>
        <dc:subject>Network data</dc:subject>
            
    <content type="html">    &lt;p&gt;This is the native binary format of DisPerSE. Functions for reading and writing &lt;em&gt;NDnet&lt;/em&gt; format  in &lt;em&gt;C&lt;/em&gt; can be found within the file &lt;code&gt;${DISPERSE_SRC}/src/C/NDnetwork.c&lt;/code&gt; (see functions &lt;em&gt;Load_NDnetwork&lt;/em&gt; and &lt;em&gt;Save_NDnetwork&lt;/em&gt;). The format may seem relatively complex, but most of it is actually optional and not used in disperse (only simplicial complexes are used in DisPerSE). To create DisPerSE input files, it is only necessary to define the highest dimensional n-simplices as a list of (n+1) vertices (see also function &lt;em&gt;CreateNetwork&lt;/em&gt;).
&lt;br /&gt;
&lt;br /&gt;
&lt;ins&gt;Note&lt;/ins&gt;: The scalar function whose MS-complex is computed by &lt;a href=&quot;http://localhost/dotclear/index.php?post/mse&quot;&gt;mse&lt;/a&gt; can be stored as an additional data field
named 'field_value' (case sensitive).
&lt;br /&gt;
&lt;ins&gt;Warning&lt;/ins&gt;: in the following, for legacy reasons, the terms &lt;em&gt;n-face&lt;/em&gt; and &lt;em&gt;n-cell&lt;/em&gt; are used indifferently to designate polygons of dimension &lt;em&gt;n&lt;/em&gt; (which are always simplexes in DisPerSE).
&lt;br /&gt;
&lt;br /&gt;
When using the C functions from Disperse, data is loaded into the following &lt;em&gt;C&lt;/em&gt; structure which is close to the actual structure of the file (see file &lt;code&gt;${DISPERSE_SRC}/src/C/NDnetwork.h&lt;/code&gt;):&lt;/p&gt;

&lt;pre&gt;typedef struct
 {
   int type; // the cell-type
   char name[255];  // name of the field
   double *data;  // value for each of the nfaces[n] n-cells
} NDnetwork_Data;&lt;/pre&gt;


&lt;pre&gt;// NDnetwork_SupData is not used in disperse ...
typedef struct
{
   int type; 
   char name[255];
   int datasize;
   char datatype[255];// a string to identity how data should be casted
   void *data;
} NDnetwork_SupData;&lt;/pre&gt;


&lt;pre&gt;typedef struct 
{
   char comment[80];
   int periodicity;
   int ndims; // the number of spatial dimensions
   int ndims_net; // number of dimension of the network itself (e.g. 2 for a sphere embedded in 3D)
   int isSimpComplex;  // 1 if network is a simplicial complex (always true in disperse)
   double *x0;  // origin of the bounding box
   double *delta;  // size of the bounding box
   int indexSize; // size of NDNET_UINT type in Bytes
   int cumIndexSize; // size of NDNET_IDCUMT type in Bytes
   char dummy[160-4*2]; // dummy data reserved for future extensions
 
   NDNET_UINT nvertex;  // total number of vertices
   float *v_coord; //vertices coodinates (X_0,Y_0,Z_0,X_1,Y_1,...,Z_nvertex-1)
   
   NDNET_UINT *nfaces; // number of cells of a given type t is given by nfaces[t]
 
   int *haveVertexFromFace; // haveVertexFromFace[n] is 1 if we have an explicit definition of the n-cells (at least one type of cell must be defined).
   NDNET_IDCUMT **f_numVertexIndexCum;// cumulative number of vertice in the t-cells, NULL when cells are simplexes (isSimpComplex=1)
   NDNET_UINT **f_vertexIndex; // list of vertices defining the n-cells is stored in f_vertexIndex[n], all vertices being enumerated for each cell (the indices of the vertices in the kth n-cell start at f_vertexIndex[n][(n+1)*k] )
   // see also macro  NUM_VERTEX_IN_FACE(net,type,face) and VERTEX_IN_FACE(net,type,face)
 
   //This may be computed internally within DisPerSE but does not need to be defined explicitely
   int *haveFaceFromVertex; // haveFaceFromVertex[n] is 1 if we have an explicit list of all the n-cells that contain each vertex (used to navigate within the network)
   NDNET_IDCUMT **v_numFaceIndexCum; // cumulative number of t-cells a vertex v belongs to
   NDNET_UINT **v_faceIndex; // indices of the t-cells in the co-boundary of v ( the list of n-cells of vertex k starts at v_faceIndex[n][net-&amp;gt;v_numFaceIndexCum[n][k]] and ends at v_faceIndex[n][net-&amp;gt;v_numFaceIndexCum[n][k+1]] )
   // see also macro  NUM_FACE_IN_VERTEX(net,type,vertex) and  FACE_IN_VERTEX(net,type,vertex)
   
   // This can become extremely memory heavy ... NOT used in DisPerSE
   int **haveFaceFromFace; // haveFaceFromFace[k][n] is 1 if we have an explicit list of all the n-cells that have a boundary/co-boundary relation with each k-cell (used to navigate within the network)
   NDNET_IDCUMT ***f_numFaceIndexCum; //  cumulative number of n-cells having a boundary / co-boundary relation with each k-cell: f_numFaceIndexCum[k][n]
   NDNET_UINT ***f_faceIndex; // indices of the cells (similar to v_faceIndex)
   // see also macro NUM_FACE_IN_FACE(net,ref_type,ref_face,type) and FACE_IN_FACE(net,ref_type,ref_face,type)
   
   int haveVFlags;  // do we have flags associated to each vertex ?
   int *haveFFlags;  // do we have flags associated to each n-cell ?
   unsigned char *v_flag; // nvertex flag values (1 for each vertex) or NULL 
   unsigned char **f_flag; // nfaces[n] flag values (1 of each n-cell) or NULL

   int ndata; // number of additional data fields.
   NDnetwork_Data *data; // array of all additionnal data (data in total)
   
   int nsupData;
   NDnetwork_SupData *supData;

} NDnetwork;&lt;/pre&gt;

&lt;p&gt;&lt;br /&gt;
&lt;br /&gt;
The &lt;em&gt;NDnet&lt;/em&gt; binary format is organized as follows (blocks are delimited by &lt;em&gt;dummy&lt;/em&gt; variables indicating the size of the blocks for FORTRAN compatibility, but they are ignored in C):
&lt;br /&gt;
&lt;br /&gt;&lt;/p&gt;
&lt;table style=&quot;margin: 1em auto 1em auto&quot;&gt;&lt;caption&gt;NDnet binary format&lt;/caption&gt;&lt;tr  valign=&quot;top&quot;&gt;&lt;td width=&quot;10%&quot;&gt;field&lt;/td&gt;&lt;td width=&quot;10%&quot;&gt;type&lt;/td&gt;&lt;td width=&quot;10%&quot;&gt;size&lt;/td&gt;&lt;td width=&quot;65%&quot;&gt;comment&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt; for FORTRAN compatibility&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;tag&lt;/td&gt;&lt;td&gt;char(1B)&lt;/td&gt;&lt;td&gt;16&lt;/td&gt;&lt;td&gt;identifies the file type. Value : &quot;NDNETWORK&quot;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;ndims&lt;/td&gt;&lt;td&gt;int(4B)&lt;/td&gt;&lt;td&gt;1&lt;/td&gt;&lt;td&gt;number of dimensions of the embedding space&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;ndims_net&lt;/td&gt;&lt;td&gt;int(4B)&lt;/td&gt;&lt;td&gt;1&lt;/td&gt;&lt;td&gt;ndims spanned by the network (=ndims by default)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;comment&lt;/td&gt;&lt;td&gt;char(1B)&lt;/td&gt;&lt;td&gt;80&lt;/td&gt;&lt;td&gt;a comment on the file (string)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt; periodicity&lt;/td&gt;&lt;td&gt;int(4B)&lt;/td&gt;&lt;td&gt;1&lt;/td&gt;&lt;td&gt; 0=non periodic, if p^th bit is set, boundary are periodic along dimension p&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;isSimpComplex&lt;/td&gt;&lt;td&gt;int(4B)&lt;/td&gt;&lt;td&gt;1&lt;/td&gt;&lt;td&gt; 1 if network is made of simplices (must be 1 for DisPerSE)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;x0&lt;/td&gt;&lt;td&gt;double(8B)&lt;/td&gt;&lt;td&gt;ndims&lt;/td&gt;&lt;td&gt;origin of bounding box&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;delta&lt;/td&gt;&lt;td&gt;double(8B)&lt;/td&gt;&lt;td&gt;ndims&lt;/td&gt;&lt;td&gt;size of bounding box&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;index_size&lt;/td&gt;&lt;td&gt;int(4B)&lt;/td&gt;&lt;td&gt;1&lt;/td&gt;&lt;td&gt;size of NDNET_UINT integer format in Bytes&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;cumindex_size&lt;/td&gt;&lt;td&gt;int(4B)&lt;/td&gt;&lt;td&gt;1&lt;/td&gt;&lt;td&gt;size of NDNET_IDCUMT integer format in Bytes&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;dummy_ext&lt;/td&gt;&lt;td&gt;char(1B)&lt;/td&gt;&lt;td&gt;152&lt;/td&gt;&lt;td&gt;dummy data reserved for future extensions&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;nvertex&lt;/td&gt;&lt;td&gt;NDNET_UINT&lt;/td&gt;&lt;td&gt;1&lt;/td&gt;&lt;td&gt;number of vertices&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;v_coords&lt;/td&gt;&lt;td&gt;float(4B)&lt;/td&gt;&lt;td&gt;ndims&amp;timesnvertex&lt;/td&gt;&lt;td&gt;coordinates of the vertices [X0,Y0, ...]&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;nfaces&lt;/td&gt;&lt;td&gt;NDNET_UINT&lt;/td&gt;&lt;td&gt;ndims+1&lt;/td&gt;&lt;td&gt;number of cells of each type (N0,N1,...)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;haveVertexFromFace&lt;/td&gt;&lt;td&gt;int(4B)&lt;/td&gt;&lt;td&gt;ndims+1&lt;/td&gt;&lt;td&gt;are n-cells explicitly defined ? (0=no, 1=yes)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;*&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;*&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;*&lt;/td&gt;&lt;td style=&quot;background:red&quot;&gt; next 3 lines are repeated for each (ndims+1) possible cells type , only if haveVertexFromFace[n] is true.&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:lightblue&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;f_vertexIndex[n]&lt;/td&gt;&lt;td&gt;NDNET_UINT&lt;/td&gt;&lt;td&gt;(n+1)&amp;timesnfaces[n]&lt;/td&gt;&lt;td style=&quot;background:lightblue&quot;&gt;list of (n+1) vertex indices for each n-cell&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:lightblue&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;haveFaceFromVertex&lt;/td&gt;&lt;td&gt;int(4B)&lt;/td&gt;&lt;td&gt;ndims+1&lt;/td&gt;&lt;td&gt;are n-cells in the co-boundary of each vertex explicitly defined ? (0=no, 1=yes)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;*&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;*&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;*&lt;/td&gt;&lt;td style=&quot;background:red&quot;&gt; next 6 lines are repeated for each (ndims+1) possible cells type, only if haveFaceFromVertex[n] is true.&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:lightblue&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;numFaceIndexCum[n]&lt;/td&gt;&lt;td&gt;NDNET_IDCUMT&lt;/td&gt;&lt;td&gt;nvertex+1&lt;/td&gt;&lt;td style=&quot;background:lightblue&quot;&gt;cumulative count of n-cells on vertices co-boundary&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:lightblue&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:lightblue&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;v_faceIndex[n]&lt;/td&gt;&lt;td&gt;NDNET_UINT&lt;/td&gt;&lt;td&gt;numFaceIndexCum[n]&lt;/td&gt;&lt;td style=&quot;background:lightblue&quot;&gt;list of n-cells on the co-boundary of each vertex (vertex i have numFaceIndexCum[n][i+1]-numFaceIndexCum[n][i] of them)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:lightblue&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;haveFaceFromFace&lt;/td&gt;&lt;td&gt;int(4B)&lt;/td&gt;&lt;td&gt;(ndims+1)^2&lt;/td&gt;&lt;td&gt;are n-cells in the co-boundary of each k-cell explicitly defined ? (0=no, 1=yes)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;*&lt;/td&gt;&lt;td&gt;*&lt;/td&gt;&lt;td&gt;*&lt;/td&gt;&lt;td style=&quot;background:red&quot;&gt; This section describes boundary relation between n-cells and k-cells. It is usually empty in DisPerSE (see NDnetwork.c for details) so SKIP IT :)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;haveVFlags&lt;/td&gt;&lt;td&gt;int(4B)&lt;/td&gt;&lt;td&gt;1&lt;/td&gt;&lt;td&gt; 1 if vertex flags are defined&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;*&lt;/td&gt;&lt;td&gt;*&lt;/td&gt;&lt;td&gt;*&lt;/td&gt;&lt;td style=&quot;background:red&quot;&gt; next 3 lines are skipped if haveVFlags=0&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:lightblue&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;v_flag&lt;/td&gt;&lt;td&gt;uchar(1B)&lt;/td&gt;&lt;td&gt;nvertex&lt;/td&gt;&lt;td style=&quot;background:lightblue&quot;&gt; value of the flags for vertices&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:lightblue&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;haveFFlags&lt;/td&gt;&lt;td&gt;int(4B)&lt;/td&gt;&lt;td&gt;ndims+1&lt;/td&gt;&lt;td&gt; 1 if flags are defined for n-cells&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;*&lt;/td&gt;&lt;td&gt;*&lt;/td&gt;&lt;td&gt;*&lt;/td&gt;&lt;td style=&quot;background:red&quot;&gt; next 3 lines are repeated for each n-cell such that haveFFlags[n]=1&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:lightblue&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;f_flag[n]&lt;/td&gt;&lt;td&gt;uchar(1B)&lt;/td&gt;&lt;td&gt;nfaces[n]&lt;/td&gt;&lt;td style=&quot;background:lightblue&quot;&gt; value of the flags for n-cells&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:lightblue&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;ndata&lt;/td&gt;&lt;td&gt;int(4B)&lt;/td&gt;&lt;td&gt;1&lt;/td&gt;&lt;td&gt; total number of additional fields&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;*&lt;/td&gt;&lt;td&gt;*&lt;/td&gt;&lt;td&gt;*&lt;/td&gt;&lt;td style=&quot;background:red&quot;&gt; next 7 lines are repeated for each additional field (&amp;timesndata)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:lightblue&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;type&lt;/td&gt;&lt;td&gt;int(4B)&lt;/td&gt;&lt;td&gt;1&lt;/td&gt;&lt;td style=&quot;background:lightblue&quot;&gt; the type of cells (0=vertex, n = n-simplex)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;name&lt;/td&gt;&lt;td&gt;char(1B)&lt;/td&gt;&lt;td&gt;255&lt;/td&gt;&lt;td style=&quot;background:lightblue&quot;&gt; name of the supplementary data&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:lightblue&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:lightblue&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;data&lt;/td&gt;&lt;td&gt;double(8B)&lt;/td&gt;&lt;td&gt;N&lt;/td&gt;&lt;td style=&quot;background:lightblue&quot;&gt; data associated to cells or vertices. N is nfaces[n] or nvertex depending on the type value.&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:lightblue&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;/table&gt;</content>
    
    

    
      </entry>
    
  <entry>
    <title>Network files</title>
    <link href="http://localhost/dotclear/index.php?post/networks" rel="alternate" type="text/html"
    title="Network files" />
    <id>urn:md5:155e7207695dd84931742192590c0813</id>
    <published>2425-05-09T06:54:00+01:00</published>
    <updated>2012-07-18T11:18:40+02:00</updated>
    <author><name>thierry sousbie</name></author>
        <dc:subject>Network data</dc:subject>
            
    <content type="html">    &lt;p&gt;This file type is designed to store any kind of unstructured network (such as delaunay tesselations or more generally cell complexes). In DisPerSE, its usage is restricted to networks of simplices though, and it is mainly used to store ascending and descending manifolds (voids and walls for instance) and persistence pairs as output by &lt;a href=&quot;http://localhost/dotclear/index.php?post/mse&quot;&gt;mse&lt;/a&gt; or Delaunay tessellations as output by &lt;a href=&quot;http://localhost/dotclear/index.php?post/Usage&quot;&gt;delaunay_nD&lt;/a&gt; (&lt;a href=&quot;http://localhost/dotclear/index.php?post/skeleton-formats&quot;&gt;skeleton files&lt;/a&gt; can also be converted to networks using &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv&quot;&gt;skelconv&lt;/a&gt; -to NDnet&lt;/em&gt;). Within network files, networks are represented by setz of vertices and cells of any dimension. A n-cell is a cell of dimension &lt;em&gt;n&lt;/em&gt;, which is therefore described, in the case of a simplicial network, as a set of &lt;em&gt;n+1&lt;/em&gt; vertices index. In the case of a simplicial &lt;em&gt;complex&lt;/em&gt;, only the highest dimensional cell have to be explicitly given, but other type of cells may also be specified. Indeed, extended manifolds for instance are not described as complexes: an ascending 0-manifold in 3D is a set of tetrahedrons (3-cells), triangles (representing ascending 1-manifolds on its boundary), segments (ascending 2-manifolds on its boundary) and vertices (ascending 3-manifolds / critical points). Note that additional information can also be associated to each type of cell (see below).
&lt;br /&gt;
The base network format is &lt;a href=&quot;http://localhost/dotclear/index.php?post/NDnet-format&quot;&gt;NDnet&lt;/a&gt; which is used internally, but this format may be converted to several other more or less complex formats of network files adapted to different applications (see option &lt;em&gt;-to&lt;/em&gt; in program &lt;a href=&quot;http://localhost/dotclear/index.php?post/netconv&quot;&gt;netconv&lt;/a&gt;, a list of available formats is displayed when running the program without argument).&lt;br /&gt;
&lt;br /&gt;
&lt;br /&gt;
&lt;strong&gt;&lt;ins&gt;Available formats&lt;/ins&gt;&lt;/strong&gt;:
&lt;br /&gt;
&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/NDnet-format&quot;&gt;NDnet&lt;/a&gt;&lt;/strong&gt; (Read / Write):&lt;br /&gt; This is the format of the network files created or red by &lt;a href=&quot;http://localhost/dotclear/index.php?post/mse&quot;&gt;mse&lt;/a&gt;. It is a relatively complex binary format (it is actually more complex than needed as it is designed to store generic non-simplicial networks) that contains all the information on the geometry and topology of unstructured networks as well as additional data associated to each type of cells.&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/NDnet_ascii-format&quot;&gt;NDnet_ascii&lt;/a&gt;&lt;/strong&gt; (Read / Write):&lt;br /&gt; This ASCII format contains the same amount of information as &lt;em&gt;NDnet&lt;/em&gt; files, but restricted to simplicial networks. It is easy to read and write so it may be used to write reasonably sized networks used as input for &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse&quot;&gt;mse&lt;/a&gt;&lt;/em&gt;.&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/ply-format&quot;&gt;PLY&lt;/a&gt;&lt;/strong&gt; and &lt;strong&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/ply-format&quot;&gt;PLY_ascii&lt;/a&gt;&lt;/strong&gt; (Read binary only / Write):&lt;br /&gt; This is a rather popular and simple binary or ASCII format that can be used as an interface with other software.&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/vtk-network-format&quot;&gt;vtk&lt;/a&gt;&lt;/strong&gt;, &lt;strong&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/vtk-network-format&quot;&gt;vtk_ascii&lt;/a&gt;&lt;/strong&gt;, &lt;strong&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/vtk-network-format&quot;&gt;vtu&lt;/a&gt;&lt;/strong&gt; and &lt;strong&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/vtk-network-formats&quot;&gt;vtu_ascii&lt;/a&gt;&lt;/strong&gt; (Write only):&lt;br /&gt; These formats are binary and ASCII legacy and XML &lt;a href=&quot;http://www.vtk.org/&quot;&gt;VTK&lt;/a&gt; formats that are readable by several 3D visualization tools, such as &lt;a href=&quot;https://wci.llnl.gov/codes/visit/&quot;&gt;VisIt&lt;/a&gt; or &lt;a href=&quot;http://www.paraview.org/&quot;&gt;ParaView&lt;/a&gt; for instance.&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;
&lt;br /&gt;
&lt;strong&gt;&lt;ins&gt;Additional data&lt;/ins&gt;&lt;/strong&gt;: In addition to the topology and geometry of the network, arbitrary additional information may be associated to each type of cell. Run &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/netconv&quot;&gt;netconv&lt;/a&gt; filename -info&lt;/em&gt; for a list of additional data available in the network file &lt;em&gt;filename&lt;/em&gt;. By default, the name of additional data added by &lt;a href=&quot;http://localhost/dotclear/index.php?post/mse&quot;&gt;mse&lt;/a&gt; is relatively explicit, it includes (see also the &lt;em&gt;additional data&lt;/em&gt; section of the &lt;a href=&quot;http://localhost/dotclear/index.php?post/skeleton-formats&quot;&gt;skeleton file format&lt;/a&gt; description):
&lt;br /&gt;
&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;field_value&lt;/strong&gt; / &lt;strong&gt;log_field_value&lt;/strong&gt; :&lt;br /&gt; The value of the field and its logarithm. The tag &lt;em&gt;field_value&lt;/em&gt; corresponds to the input function for &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse&quot;&gt;mse&lt;/a&gt;&lt;/em&gt;, whose Morse-Smale complex is to be computed.&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;cell&lt;/strong&gt;:&lt;br /&gt; The &lt;em&gt;type&lt;/em&gt; and &lt;em&gt;index&lt;/em&gt; of a cell in the original network (prefix may be added). The value is a double precision floating number whose integer part is the index of the cell and decimal part its type. For instance, the 156th vertex (i.e. 0-cell) in the cell complex is represented as 156.0, while the 123th tetrahedron is 123.3. Note that the index of the 0-cell correpond to the index of the pixel / vertices in the original network from which the skeleton was computed.&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;type&lt;/strong&gt;:&lt;br /&gt; This usually corresponds to the critical index of a critical point (for instance, vertices of persistence pairs networks), or the type of a persistence pair (i.e. the minimum critical index of the CP in the pair, for segments of persistence pairs networks).&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;index&lt;/strong&gt;&lt;br /&gt; Usually the index of a vertex (e.g. for persistence pairs, additional segment data tagged &lt;em&gt;up_index&lt;/em&gt; and &lt;em&gt;down_index&lt;/em&gt; correspond to the indices of the vertices with lowest and highest critical index in the persistence pair respectively).&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;persistence&lt;/strong&gt; / &lt;strong&gt;persistence_ratio&lt;/strong&gt; / &lt;strong&gt;persistence_nsigmas&lt;/strong&gt; :&lt;br /&gt; The persistence (expressed as a difference, ratio or in &lt;em&gt;number of sigmas&lt;/em&gt;) of the persistence pair containing the corresponding critical point. A negative or null value indicates that &lt;em&gt;persistence&lt;/em&gt; is not relevant to this particular cell.&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;parent&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;parent_index&lt;/strong&gt; / &lt;strong&gt;parent_log_index&lt;/strong&gt; (vertices only):&lt;br /&gt; For persistence pairs type networks, for each vertex representing an extremum  (i.e. minima and maxima), the index of the vertex that corresponds to the other extremum into which it would be merged if its persistence pair was canceled (indices start at 0). This can be used to reconstruct the tree of the hierarchy of maxima and minima. The value is -1 for non extrema critical points. The difference between the two versions is that the second (&lt;em&gt;parent_log_index&lt;/em&gt;) is the hierarchy computed from the logarithm of the field. The second version is useful only for discrete point samples whose MS-complex is obtained from the delaunay tessellation computed with &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/Usage&quot;&gt;delaunay_nD&lt;/a&gt;&lt;/em&gt;.  Practically, &lt;em&gt;parent_log_index&lt;/em&gt; can be used whenever persistence pairs are cancelled in order of increasing ratio (option &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#nsig&quot;&gt;-nsig&lt;/a&gt;&lt;/em&gt; in &lt;em&gt;mse&lt;/em&gt;), and &lt;em&gt;parent_index&lt;/em&gt; whenever persistence pairs are cancelled in order of increasing difference  (option &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#cut&quot;&gt;-cut&lt;/a&gt;&lt;/em&gt; in &lt;em&gt;mse&lt;/em&gt;).&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;a name=&quot;source&quot;&gt;&lt;/a&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;source_cell&lt;/strong&gt; / &lt;strong&gt;source_index&lt;/strong&gt;:&lt;br /&gt; For networks representing manifolds (voids, walls, ... obtained with option &lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#dumpmanifolds&quot;&gt;-dumpManifolds&lt;/a&gt; of mse), this represents for each simplex the critical point from which the manifold it belongs to originates (for instance, the minimum corresponding to a void, or the saddle point corresponding to a filament). In &lt;em&gt;source_cell&lt;/em&gt;, the critical points is represented by its cell in the initial cell complex (see &lt;strong&gt;cell&lt;/strong&gt; above), while &lt;em&gt;source_index&lt;/em&gt; gives the index of the critical point in the &lt;a href=&quot;http://localhost/dotclear/index.php?post/skeleton-formats&quot;&gt;skeleton&lt;/a&gt; file or persistence pair network obtained with &lt;a href=&quot;http://localhost/dotclear/index.php?post/mse&quot;&gt;mse&lt;/a&gt; (options &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#dumparcs&quot;&gt;-dumpArcs&lt;/a&gt;&lt;/em&gt; and &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse#ppairs&quot;&gt;-ppairs&lt;/a&gt;&lt;/em&gt;). See also &lt;a href=&quot;http://localhost/dotclear/index.php?post/Example-1#voids&quot;&gt;here&lt;/a&gt; and &lt;a href=&quot;http://localhost/dotclear/index.php?post/Example-2#walls&quot;&gt;there&lt;/a&gt; in the tutorial section.&lt;/li&gt;
&lt;/ul&gt;</content>
    
    

    
      </entry>
  
</feed>