<?xml version="1.0" encoding="utf-8"?><feed xmlns="http://www.w3.org/2005/Atom"
  xmlns:dc="http://purl.org/dc/elements/1.1/"
  xmlns:wfw="http://wellformedweb.org/CommentAPI/"
  xml:lang="en">
  
  <title type="html">DisPerSE - persistent structures identification - Field data</title>
  <subtitle type="html">Automatic identification of persistent structures in 2D or 3D.
keywords: Morse complex, topology, peak, void, source, wall, filament, cosmic web, cosmology, structure identification.</subtitle>
  <link href="http://localhost/dotclear/index.php?feed/category/Field-I-O/atom" rel="self" type="application/atom+xml"/>
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
    <title>vtk field format</title>
    <link href="http://localhost/dotclear/index.php?post/vtk-field-format" rel="alternate" type="text/html"
    title="vtk field format" />
    <id>urn:md5:815da7c1d79dc8027044ce918cca30b3</id>
    <published>2430-01-06T06:56:00+01:00</published>
    <updated>2012-05-28T16:33:34+02:00</updated>
    <author><name>thierry sousbie</name></author>
        <dc:subject>Field data</dc:subject>
            
    <content type="html">    &lt;p&gt;VTK formats are developed for the Visualization Tool Kit library (&lt;a href=&quot;http://localhost/dotclear/index.php?post/vtk-formats&quot;&gt;VTK&lt;/a&gt;) and can be used for 3D visualization with software such as &lt;a href=&quot;https://wci.llnl.gov/codes/visit/&quot;&gt;VisIt&lt;/a&gt; or &lt;a href=&quot;http://www.paraview.org/&quot;&gt;ParaView&lt;/a&gt;. Only regular grids field types can be converted to VTK (for particle type fields, convert to &lt;em&gt;NDnet&lt;/em&gt; &lt;a href=&quot;http://localhost/dotclear/index.php?post/networks&quot;&gt;unstructured network&lt;/a&gt; format first and then use &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/netconv&quot;&gt;netconv&lt;/a&gt; -to vtk&lt;/em&gt;). Regular grids are stored as &lt;em&gt;VTK image data&lt;/em&gt; and can be output in fours different VTK formats:
&lt;br /&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;vtk&lt;/strong&gt;: the legacy format&lt;/li&gt;
&lt;li&gt;-&lt;strong&gt;vtk_ascii&lt;/strong&gt;: ASCII version of the vtk format&lt;/li&gt;
&lt;li&gt;-&lt;strong&gt;vti&lt;/strong&gt;: a more recently developed XML version of the vtk format,&lt;/li&gt;
&lt;li&gt;-&lt;strong&gt;vti_ascii&lt;/strong&gt;: ASCII version of the vtu format&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;
&lt;br /&gt;
The specifications for these formats can be found in this &lt;a href=&quot;http://localhost/dotclear/index.php?post/www.vtk.org/VTK/img/file-formats.pdf&quot; title=&quot;VTK formats PDF&quot;&gt;PDF&lt;/a&gt; file. See also &lt;a href=&quot;http://www.cacr.caltech.edu/~slombey/asci/vtk/vtk_formats.simple.html&quot;&gt;here&lt;/a&gt; and &lt;a href=&quot;http://mathema.tician.de/node/430&quot; title=&quot;there&quot;&gt;there&lt;/a&gt; for additional information.
&lt;br /&gt;
&lt;br /&gt;
&lt;br /&gt;
&lt;br /&gt;&lt;/p&gt;</content>
    
    

    
      </entry>
    
  <entry>
    <title>SDL-image format</title>
    <link href="http://localhost/dotclear/index.php?post/SDL-image" rel="alternate" type="text/html"
    title="SDL-image format" />
    <id>urn:md5:cdd6677cf033684304583dc51bad63db</id>
    <published>2430-01-05T06:56:00+01:00</published>
    <updated>2012-05-28T16:27:56+02:00</updated>
    <author><name>thierry sousbie</name></author>
        <dc:subject>Field data</dc:subject>
            
    <content type="html">    &lt;p&gt;This input format is used to read 2D regular grids encoded in a popular picture format readable by &lt;a href=&quot;http://www.libsdl.org/projects/SDL_image/docs/SDL_image.html&quot; title=&quot;SDL-image&quot;&gt;SDL-image&lt;/a&gt; library (jpg, png, gif, ...). When using this format as input to &lt;a href=&quot;http://localhost/dotclear/index.php?post/mse&quot;&gt;mse&lt;/a&gt;, the field whose topology is computed is the luminosity L of the image, a weighted average of the RED, GREEN and BLUE components of the image:
&lt;br /&gt;
&lt;br /&gt;&lt;/p&gt;

&lt;pre&gt;L = 0.2989*RED + 0.5870*GREEN + 0.1140*BLUE&lt;/pre&gt;

&lt;p&gt;&lt;br /&gt;
&lt;br /&gt;
&lt;br /&gt;
&lt;br /&gt;
&lt;br /&gt;
&lt;br /&gt;
&lt;br /&gt;
&lt;br /&gt;
&lt;br /&gt;
&lt;br /&gt;&lt;/p&gt;</content>
    
    

    
      </entry>
    
  <entry>
    <title>survey_ascii format</title>
    <link href="http://localhost/dotclear/index.php?post/survey_ascii-format" rel="alternate" type="text/html"
    title="survey_ascii format" />
    <id>urn:md5:0642003343a6eb7fa09f2c323576fd94</id>
    <published>2430-01-04T06:56:00+01:00</published>
    <updated>2012-05-28T17:16:07+02:00</updated>
    <author><name>thierry sousbie</name></author>
        <dc:subject>Field data</dc:subject>
            
    <content type="html">    &lt;p&gt;This is a simple ASCII format whose main purpose is to easily read coordinates of discretely sampled astrophysical surveys (e.g. such as an &lt;a href=&quot;http://www.sdss.org/&quot; title=&quot;SDSS&quot;&gt;SDSS&lt;/a&gt; galaxy catalog). It should be considered when using spherical coordinates systems or if distance is measured by redshift for instance.
&lt;br /&gt;
&lt;br /&gt;
In this format, particles properties are encoded in an ASCII array where each row corresponds to one particle and each column to one property. Each column must have a name defined in the header (the first line of the file, starting with character &lt;em&gt;#&lt;/em&gt;). An survey_ascii file may look like this:
&lt;br /&gt;
&lt;br /&gt;&lt;/p&gt;

&lt;pre&gt;# ra dec z my_field
+2.115401e+02	+6.313433e-01	+2.666800e-01  1	
+2.115633e+02	+7.550188e-01	+1.259900e-01  0	
+2.115687e+02	+8.108763e-01	+3.646600e-01  0	
+2.117158e+02	+6.393598e-01	+1.143600e-01  1	
+2.116826e+02	+6.528485e-01	+2.455700e-01  1	
+2.116993e+02	+6.509297e-01	+1.199000e-01  0	
+2.115738e+02	+7.772653e-01	+3.240600e-01  0	
+2.116198e+02	+6.950604e-01	+1.987300e-01  0	
+2.116773e+02	+7.085776e-01	+2.561900e-01  1&lt;/pre&gt;

&lt;p&gt;&lt;br /&gt;
&lt;br /&gt;
The name of a column defines its role if it matches one of the following keywords:
&lt;br /&gt;
&lt;br /&gt;&lt;/p&gt;
&lt;table border=&quot;1&quot;&gt;&lt;caption&gt; header keywords&lt;/caption&gt;&lt;tr align=&quot;center&quot;&gt;&lt;td width=&quot;20%&quot; style=&quot;background:grey&quot;&gt;Keyword&lt;/td&gt;&lt;td width=&quot;40%&quot; style=&quot;background:grey&quot;&gt;Meaning&lt;/td&gt;&lt;/tr&gt;
&lt;tr align=&quot;center&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt; px &lt;/td&gt;&lt;td&gt; X coordinate of the particle&lt;/td&gt;&lt;/tr&gt;
&lt;tr align=&quot;center&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt; py &lt;/td&gt;&lt;td&gt; Y coordinate of the particle&lt;/td&gt;&lt;/tr&gt;
&lt;tr align=&quot;center&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt; pz &lt;/td&gt;&lt;td&gt; Z coordinate of the particle&lt;/td&gt;&lt;/tr&gt;
&lt;tr align=&quot;center&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt; vx &lt;/td&gt;&lt;td&gt; X component of the velocity&lt;/td&gt;&lt;/tr&gt;
&lt;tr align=&quot;center&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt; vy &lt;/td&gt;&lt;td&gt; Y component of the velocity&lt;/td&gt;&lt;/tr&gt;
&lt;tr align=&quot;center&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt; vz &lt;/td&gt;&lt;td&gt; Z component of the velocity&lt;/td&gt;&lt;/tr&gt;
&lt;tr align=&quot;center&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt; id &lt;/td&gt;&lt;td&gt; an index associated to the particle&lt;/td&gt;&lt;/tr&gt;
&lt;tr align=&quot;center&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt; ra &lt;/td&gt;&lt;td&gt; The &lt;a href=&quot;http://en.wikipedia.org/wiki/Right_ascension&quot;&gt;right ascension&lt;/a&gt; of the particle&lt;/td&gt;&lt;/tr&gt;
&lt;tr align=&quot;center&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt; dec &lt;/td&gt;&lt;td&gt; The &lt;a href=&quot;http://en.wikipedia.org/wiki/Declination&quot;&gt;declination&lt;/a&gt; of the particle&lt;/td&gt;&lt;/tr&gt;
&lt;tr align=&quot;center&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt; z &lt;/td&gt;&lt;td&gt; The &lt;a href=&quot;http://en.wikipedia.org/wiki/Redshift&quot;&gt;redshift&lt;/a&gt; of the particle.&lt;/td&gt;&lt;/tr&gt;
&lt;/table&gt;

&lt;p&gt;&lt;br /&gt;
&lt;br /&gt;
When particles coordinates are defined with &lt;em&gt;ra&lt;/em&gt;, &lt;em&gt;dec&lt;/em&gt; or &lt;em&gt;z&lt;/em&gt; (redshift), DisPerSE automatically transform them into cartesian coordinates using the standard LCDM model to compute distances (Omega_m=0.27, Omega_L=0.73 ). This transformation can be inverted on the output skeletons and networks using options &lt;em&gt;-toRaDecZ&lt;/em&gt; and &lt;em&gt;-toRaDecDist&lt;/em&gt; of &lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv&quot;&gt;skelconv&lt;/a&gt; and &lt;a href=&quot;http://localhost/dotclear/index.php?post/netconv&quot;&gt;netconv&lt;/a&gt;.&lt;/p&gt;</content>
    
    

    
      </entry>
    
  <entry>
    <title>FITS format</title>
    <link href="http://localhost/dotclear/index.php?post/FITS-format" rel="alternate" type="text/html"
    title="FITS format" />
    <id>urn:md5:4a123cf8572bc22997acb2ed219cd137</id>
    <published>2430-01-03T06:56:00+01:00</published>
    <updated>2012-05-28T16:21:06+02:00</updated>
    <author><name>thierry sousbie</name></author>
        <dc:subject>Field data</dc:subject>
            
    <content type="html">    &lt;p&gt;This field format is used to read 2D or 3D regular grids of arbitrary dimensions in the popular Flexible Image Transport System (&lt;a href=&quot;http://en.wikipedia.org/wiki/FITS&quot; title=&quot;FITSS format&quot;&gt;FITS&lt;/a&gt;) image format. A popular library for reading and writing FITS image in &lt;em&gt;C&lt;/em&gt; or &lt;em&gt;FORTRAN&lt;/em&gt; is &lt;a href=&quot;http://heasarc.gsfc.nasa.gov/docs/software/fitsio/fitsio.html&quot; title=&quot;FITS IO&quot;&gt;CFITSIO&lt;/a&gt;, which must be installed for DisPerSE to be able to use this kind of files.
&lt;br /&gt;
&lt;br /&gt;
This file format is also used for storing &lt;a href=&quot;http://healpix.jpl.nasa.gov/&quot; title=&quot;HEALPIX&quot;&gt;HEALPIX&lt;/a&gt; tesselations of the sphere. DisPerSE will automatically detect if the FITS file contains a HEALPIX tesselation.
&lt;br /&gt;
&lt;br /&gt;
&lt;br /&gt;
&lt;br /&gt;
&lt;br /&gt;
&lt;br /&gt;
&lt;br /&gt;
&lt;br /&gt;
&lt;br /&gt;
&lt;br /&gt;
&lt;br /&gt;&lt;/p&gt;</content>
    
    

    
      </entry>
    
  <entry>
    <title>NDfield format</title>
    <link href="http://localhost/dotclear/index.php?post/NDfield-format" rel="alternate" type="text/html"
    title="NDfield format" />
    <id>urn:md5:3aaeb963f82937bce4a9fd6addb5f958</id>
    <published>2430-01-02T06:56:00+01:00</published>
    <updated>2012-05-28T15:56:08+02:00</updated>
    <author><name>thierry sousbie</name></author>
        <dc:subject>Field data</dc:subject>
            
    <content type="html">    &lt;p&gt;This is the native binary format of DisPerSE. Functions for reading and writing &lt;em&gt;NDfield&lt;/em&gt; format  in &lt;em&gt;C&lt;/em&gt; can be found within the file &lt;code&gt;${DISPERSE_SRC}/src/C/NDfield.c&lt;/code&gt; (see functions &lt;em&gt;Load_NDfield&lt;/em&gt; and &lt;em&gt;Save_NDfield&lt;/em&gt;).
&lt;br /&gt;
&lt;br /&gt;
When using the C functions from Disperse, data is loaded into the following &lt;em&gt;C&lt;/em&gt; structure which is close to the actual structure of the file (see file &lt;code&gt;${DISPERSE_SRC}/src/C/NDfield.h&lt;/code&gt;):&lt;/p&gt;


&lt;pre&gt;#define ND_CHAR   (1&amp;lt;&amp;lt;0)
#define ND_UCHAR  (1&amp;lt;&amp;lt;1)
#define ND_SHORT  (1&amp;lt;&amp;lt;2)
#define ND_USHORT (1&amp;lt;&amp;lt;3)
#define ND_INT    (1&amp;lt;&amp;lt;4)
#define ND_UINT   (1&amp;lt;&amp;lt;5)
#define ND_LONG   (1&amp;lt;&amp;lt;6)
#define ND_ULONG  (1&amp;lt;&amp;lt;7)
#define ND_FLOAT  (1&amp;lt;&amp;lt;8)
#define ND_DOUBLE (1&amp;lt;&amp;lt;9)&lt;/pre&gt;


&lt;pre&gt;typedef struct NDfield_str
{
 char comment[80];  // a comment on the data
 int dims[NDFIELD_MAX_DIMS];  // dimensions of the grid, must be [ndims,nparticles] when data represents sample particles coordinates (i.e. when fdims_index!=0)
 int ndims;  // number of dimensions of the space
 int n_dims;  // number of meaningfull values in dims array
 int fdims_index; // if 0, the field is a regular grid of dimensions dims, else the file contains the dims[0] coordinates of dims[1] particles.
 int datatype;  // type of the data (one of the ND_... defined above)
 double x0[NDFIELD_MAX_DIMS];  // origin of the bounding box
 double delta[NDFIELD_MAX_DIMS];  // extent of the bounding box
 char dummy[160];  // dummy data, for future extensions or for storing anything you want.

 void *val;  // pointer to data

 long nval;  // total number of particles (fdims_index==1) or pixels (fdims_index==0)
 int datasize;  // size in bytes of datatype type.
} NDfield;&lt;/pre&gt;


&lt;p&gt;&lt;br /&gt;
&lt;br /&gt;
The &lt;em&gt;NDfield&lt;/em&gt; binary format is organized as follows (blocks are delimited by &lt;em&gt;dummy&lt;/em&gt; variables indicating the size of the blocks for FORTRAN compatibility, but they are ignored in C):
&lt;br /&gt;
&lt;br /&gt;&lt;/p&gt;
&lt;table style=&quot;margin: 1em auto 1em auto&quot;&gt;&lt;caption&gt;NDnet binary format&lt;/caption&gt;&lt;tr &gt;&lt;td width=&quot;10%&quot; &gt;field&lt;/td&gt;&lt;td width=&quot;10%&quot; &gt;type&lt;/td&gt;&lt;td width=&quot;10%&quot;&gt;size&lt;/td&gt;&lt;td width=&quot;65%&quot; &gt;comment&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt; for FORTRAN compatibility&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;tag&lt;/td&gt;&lt;td&gt;char(1B)&lt;/td&gt;&lt;td&gt;16&lt;/td&gt;&lt;td&gt;identifies the file type. Value : &quot;NDFIELD&quot;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;ndims&lt;/td&gt;&lt;td&gt;int(4B)&lt;/td&gt;&lt;td&gt;1&lt;/td&gt;&lt;td&gt;number of dimensions of the embedding space&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;dims&lt;/td&gt;&lt;td&gt;int(4B)&lt;/td&gt;&lt;td&gt;20&lt;/td&gt;&lt;td&gt;size of the grid in pixels along each dimension, or [ndims,nparticles] if data represents particle coordinates (i.e. fdims_index=1)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;fdims_index&lt;/td&gt;&lt;td&gt;int(4B)&lt;/td&gt;&lt;td&gt;1&lt;/td&gt;&lt;td&gt; 0 if data represents a regular grid, 1 if it represents coordinates of tracer particles&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;datatype&lt;/td&gt;&lt;td&gt;int(4B)&lt;/td&gt;&lt;td&gt;1&lt;/td&gt;&lt;td&gt;type of data stored (see below)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;x0&lt;/td&gt;&lt;td&gt;double(8B)&lt;/td&gt;&lt;td&gt;20&lt;/td&gt;&lt;td&gt;origin of bounding box (first ndims val. are meaningfull)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;delta&lt;/td&gt;&lt;td&gt;double(8B)&lt;/td&gt;&lt;td&gt;20&lt;/td&gt;&lt;td&gt;size of bounding box (first ndims val. are meaningfull)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;dummy_ext&lt;/td&gt;&lt;td&gt;char(1B)&lt;/td&gt;&lt;td&gt;160&lt;/td&gt;&lt;td&gt;dummy data reserved for future extensions&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td&gt;data&lt;/td&gt;&lt;td&gt;size of datatype&lt;/td&gt;&lt;td&gt;N&lt;/td&gt;&lt;td&gt;data itself (N may be the number of pixels or ndism times the number of particles)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt;dummy&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;int(4B) &lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;1&lt;/td&gt;&lt;td style=&quot;background:grey&quot;&gt;&lt;/td&gt;&lt;/tr&gt;
&lt;/table&gt;

&lt;p&gt;&lt;br /&gt;
&lt;br /&gt;&lt;/p&gt;
&lt;table border=&quot;1&quot;&gt;&lt;caption&gt; Possible data types&lt;/caption&gt;&lt;tr align=&quot;center&quot;&gt;&lt;td width=&quot;20%&quot; style=&quot;background:grey&quot;&gt;name&lt;/td&gt;&lt;td width=&quot;20%&quot; style=&quot;background:grey&quot;&gt;size (Bytes)&lt;/td&gt;&lt;td width=&quot;20%&quot; style=&quot;background:grey&quot;&gt;type&lt;/td&gt;&lt;td width=&quot;20%&quot; style=&quot;background:grey&quot;&gt;value&lt;/td&gt;&lt;/tr&gt;
&lt;tr align=&quot;center&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt; ND_CHAR &lt;/td&gt;&lt;td&gt; 1 &lt;/td&gt;&lt;td&gt; integer &lt;/td&gt;&lt;td&gt;  1 (=1&lt;&lt;0)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  align=&quot;center&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt; ND_UCHAR &lt;/td&gt;&lt;td&gt; 1 &lt;/td&gt;&lt;td&gt; integer &lt;/td&gt;&lt;td&gt;  2 (=1&lt;&lt;1)&lt;/td&gt;&lt;/tr&gt;
&lt;tr align=&quot;center&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt; ND_SHORT &lt;/td&gt;&lt;td&gt; 2 &lt;/td&gt;&lt;td&gt; integer &lt;/td&gt;&lt;td&gt;  4 (=1&lt;&lt;2)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  align=&quot;center&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt; ND_USHORT &lt;/td&gt;&lt;td&gt; 2 &lt;/td&gt;&lt;td&gt; integer &lt;/td&gt;&lt;td&gt;  8 (=1&lt;&lt;3)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  align=&quot;center&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt; ND_INT &lt;/td&gt;&lt;td&gt; 4 &lt;/td&gt;&lt;td&gt; integer &lt;/td&gt;&lt;td&gt;  16 (=1&lt;&lt;4)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  align=&quot;center&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt; ND_UINT  &lt;/td&gt;&lt;td&gt; 4 &lt;/td&gt;&lt;td&gt; integer &lt;/td&gt;&lt;td&gt;  32 (=1&lt;&lt;5)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  align=&quot;center&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt; ND_LONG &lt;/td&gt;&lt;td&gt; 8 &lt;/td&gt;&lt;td&gt; integer &lt;/td&gt;&lt;td&gt;  64 (=1&lt;&lt;6)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  align=&quot;center&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt; ND_ULONG &lt;/td&gt;&lt;td&gt; 8 &lt;/td&gt;&lt;td&gt; integer &lt;/td&gt;&lt;td&gt;  128 (=1&lt;&lt;7)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  align=&quot;center&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt; ND_FLOAT &lt;/td&gt;&lt;td&gt; 4 &lt;/td&gt;&lt;td&gt; float &lt;/td&gt;&lt;td&gt;  256 (=1&lt;&lt;8)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  align=&quot;center&quot;&gt;&lt;td style=&quot;background:grey&quot;&gt; ND_DOUBLE &lt;/td&gt;&lt;td&gt; 8 &lt;/td&gt;&lt;td&gt; float &lt;/td&gt;&lt;td&gt;  512 (=1&lt;&lt;9)&lt;/td&gt;&lt;/tr&gt;
&lt;/table&gt;</content>
    
    

    
      </entry>
    
  <entry>
    <title>NDfield_ascii format</title>
    <link href="http://localhost/dotclear/index.php?post/NDfield_ascii-format" rel="alternate" type="text/html"
    title="NDfield_ascii format" />
    <id>urn:md5:13c96188942b97e6721220412233a137</id>
    <published>2430-01-02T06:56:00+01:00</published>
    <updated>2013-01-24T04:58:22+01:00</updated>
    <author><name>thierry sousbie</name></author>
        <dc:subject>Field data</dc:subject>
            
    <content type="html">    &lt;p&gt;A simple ASCII version of the &lt;a href=&quot;http://localhost/dotclear/index.php?post/NDfield-format&quot;&gt;NDfield&lt;/a&gt; format used to represent either uniform grids or particles coordinates. Data type is always interpreted as double precision floating point.
&lt;br /&gt;
&lt;br /&gt;&lt;/p&gt;
&lt;table border=&quot;1&quot;&gt;&lt;caption&gt; NDfield_ascii format&lt;/caption&gt;&lt;tr &gt;&lt;td  width=&quot;50%&quot; &gt;&lt;/td&gt;&lt;td  width=&quot;50%&quot; &gt;&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt;ANDFIELD	COORDS&lt;/th&gt;&lt;td&gt;header (COORDS keyword is optional, if present, values are interpreted as coordinates, or else, as pixel values)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt;[N_1 ... N_ndims]&lt;/th&gt;&lt;td&gt;size of the grid along each of the ndims dimensions (in pixels) or number of dimensions and number of particles if COORDS keyword is in the header (= [ndims,nparticles])&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt;BBOX [X0_1 ... X0_ndims] [delta_1 ... delta_ndims]&lt;/th&gt;&lt;td&gt;OPTIONAL: bounding box definition 	(origin and extent)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt;V_1 V_2&lt;/th&gt;&lt;td&gt; values (any number of values can be on the same line)&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt;V_3&lt;/th&gt;&lt;td&gt; no particular formating is required&lt;/td&gt;&lt;/tr&gt;
&lt;tr  valign=&quot;top&quot;&gt;&lt;th&gt;...&lt;/th&gt;&lt;td&gt; Repeat enough time to define a number of values equal to the number of pixels or the number of particles (if COORDS keyword is in the header)&lt;/td&gt;&lt;/tr&gt;
&lt;/table&gt;</content>
    
    

    
      </entry>
    
  <entry>
    <title>Field files</title>
    <link href="http://localhost/dotclear/index.php?post/field-formats" rel="alternate" type="text/html"
    title="Field files" />
    <id>urn:md5:c641a01169a789e3018403e8945394c8</id>
    <published>2430-01-01T06:56:00+01:00</published>
    <updated>2012-07-13T14:51:29+02:00</updated>
    <author><name>thierry sousbie</name></author>
        <dc:subject>Field data</dc:subject>
            
    <content type="html">    &lt;p&gt;This file type is designed to store scalar fields sampled over regular grids (i.e. regular images in 2D, 3D and more) or, by extension, coordinates of tracer particles sampling an underlying density field. The main field format is &lt;em&gt;NDfield&lt;/em&gt;, which allows to store any type of data (integer, simple or double precision floating point, ...) over arbitrary sized regular grids or as coordinates of tracer particles.
&lt;br /&gt;
&lt;br /&gt;
&lt;strong&gt;&lt;ins&gt;Available formats&lt;/ins&gt;&lt;/strong&gt;:
&lt;br /&gt;
&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/NDfield-format&quot;&gt;NDfield&lt;/a&gt;&lt;/strong&gt; (Read / Write):&lt;br /&gt;This is the format used internally in DisPerSE, it is efficient and generic as it can store regular grids or sample particles coordinates indifferently.&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/NDfield_ascii-format&quot;&gt;NDfield_ascii&lt;/a&gt;&lt;/strong&gt; (Read / Write):&lt;br /&gt;A simple ASCII version of the NDfield format, it is as versatile but data is always considered to be double precision floating point.&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/FITS-format&quot;&gt;FITS&lt;/a&gt;&lt;/strong&gt; (Read only):&lt;br /&gt;Used to read 2D or 3D regular grids of arbitrary dimensions in the popular Flexible Image Transport System (&lt;a href=&quot;http://en.wikipedia.org/wiki/FITS&quot; title=&quot;FITSS format&quot;&gt;FITS&lt;/a&gt;) image format. Also used for &lt;a href=&quot;http://healpix.jpl.nasa.gov/&quot; title=&quot;HEALPIX&quot;&gt;HEALPIX&lt;/a&gt; tessellations of the sphere.&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/survey_ascii-format&quot;&gt;survey_ascii&lt;/a&gt;&lt;/strong&gt;(Read only):&lt;br /&gt;A simple ASCII format whose main purpose is to easily read coordinates of discretely sampled astrophysical surveys (e.g. such as an &lt;a href=&quot;http://www.sdss.org/&quot; title=&quot;SDSS&quot;&gt;SDSS&lt;/a&gt; galaxy catalog). Try this if you use spherical coordinates or if distance is measured by redshift for instance.&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/SDL-image&quot;&gt;SDL-image&lt;/a&gt;&lt;/strong&gt; (Read only):&lt;br /&gt;This input format is used to read 2D regular grids encoded in a popular picture format readable by &lt;a href=&quot;http://www.libsdl.org/projects/SDL_image/docs/SDL_image.html&quot; title=&quot;SDL-image&quot;&gt;SDL-image&lt;/a&gt; library (jpg, png, gif, ...).&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/vtk-field-format&quot;&gt;vtk&lt;/a&gt;&lt;/strong&gt;, &lt;strong&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/vtk-field-format&quot;&gt;vtk_ascii&lt;/a&gt;&lt;/strong&gt;, &lt;strong&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/vtk-field-format&quot;&gt;vti&lt;/a&gt;&lt;/strong&gt; and &lt;strong&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/vtk-field-format&quot;&gt;vti_ascii&lt;/a&gt;&lt;/strong&gt; (Write only):&lt;br /&gt;These formats are binary and ASCII legacy and XML &lt;a href=&quot;http://www.vtk.org/&quot;&gt;VTK&lt;/a&gt; formats that are readable by several 3D visualization tools, such as &lt;a href=&quot;https://wci.llnl.gov/codes/visit/&quot;&gt;VisIt&lt;/a&gt; or &lt;a href=&quot;http://www.paraview.org/&quot;&gt;ParaView&lt;/a&gt; for instance.&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;&lt;/p&gt;
&lt;ul&gt;
&lt;li&gt;-&lt;strong&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/networks&quot;&gt;NDnet&lt;/a&gt;&lt;/strong&gt; (Write only):&lt;br /&gt; Using this format, field files representing tracer particles coordinates can be converted to &lt;a href=&quot;http://localhost/dotclear/index.php?post/networks&quot;&gt;unstructured networks&lt;/a&gt; (use &lt;em&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/netconv&quot;&gt;netconv&lt;/a&gt;&lt;/em&gt; to convert &lt;em&gt;NDnet&lt;/em&gt; files to other network formats). Regular grids conversion is not implemented yet due to the huge size of the output network.&lt;/li&gt;
&lt;/ul&gt;
&lt;p&gt;&lt;br /&gt;
&lt;br /&gt;&lt;/p&gt;</content>
    
    

    
      </entry>
  
</feed>