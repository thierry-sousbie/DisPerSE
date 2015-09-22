<?xml version="1.0" encoding="utf-8"?><feed xmlns="http://www.w3.org/2005/Atom"
  xmlns:dc="http://purl.org/dc/elements/1.1/"
  xmlns:wfw="http://wellformedweb.org/CommentAPI/"
  xml:lang="en">
  
  <title type="html">DisPerSE - persistent structures identification - Install</title>
  <subtitle type="html">Automatic identification of persistent structures in 2D or 3D.
keywords: Morse complex, topology, peak, void, source, wall, filament, cosmic web, cosmology, structure identification.</subtitle>
  <link href="http://localhost/dotclear/index.php?feed/category/Install/atom" rel="self" type="application/atom+xml"/>
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
    <title>Compiling</title>
    <link href="http://localhost/dotclear/index.php?post/Compiling" rel="alternate" type="text/html"
    title="Compiling" />
    <id>urn:md5:a9530f8f886c454715738918e902506a</id>
    <published>2102-05-08T13:52:00+01:00</published>
    <updated>2012-11-15T14:13:41+01:00</updated>
    <author><name>thierry sousbie</name></author>
        <dc:subject>Install</dc:subject>
            
    <content type="html">    Go to the &lt;q&gt;${DISPERSE_SRC}/build&lt;/q&gt; directory and type:&lt;div&gt;&amp;nbsp;
&lt;div style=&quot;text-align: center; border: 1px dotted gray;&quot;&gt;cmake ../&lt;/div&gt;
&lt;br /&gt;This will check the configuration and generate the Makefile. Read the output to know which library were found, which were not, and how to specify their path (option &lt;q&gt;-DLIBNAME_DIR=path/to/library/&lt;/q&gt; where {LIBNAME} may be QT,GSL,SDL,MATHGL or CGAL).&lt;/div&gt;&lt;div&gt;The default install path is &quot;${DISPERSE_SRC}&quot;, but it can be changed using option &lt;q&gt;-DCMAKE_INSTALL_PREFIX=PATH/TO/INSTALL/DIR&lt;/q&gt;.&lt;br /&gt;Then just compile and install with:&lt;/div&gt;&lt;div&gt;&lt;br /&gt;&lt;div style=&quot;text-align: center; border: 1px dotted gray;&quot;&gt;make install -j N&lt;/div&gt;&lt;br /&gt;where N is the number of processors to be used for compilation. One to seven executables should be created in &lt;q&gt;${DISPERSE_SRC}/bin&lt;/q&gt; (or &lt;q&gt;${CMAKE_INSTALL_PREFIX}/bin&lt;/q&gt; if set), depending on which libraries were found:&amp;nbsp;&lt;div&gt;&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/mse&quot;&gt;mse&lt;/a&gt;&lt;/li&gt;
&lt;li&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/Usage&quot;&gt;delaunay_2D&lt;/a&gt; (Optionnal)&lt;/li&gt;
&lt;li&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/Usage&quot;&gt;delaunay_3D&lt;/a&gt; (Optionnal)&lt;/li&gt;
&lt;li&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/pdview&quot;&gt;pdview&lt;/a&gt; (Optionnal)&lt;/li&gt;
&lt;li&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/skelconv&quot;&gt;skelconv&lt;/a&gt;&lt;/li&gt;
&lt;li&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/netconv&quot;&gt;netconv&lt;/a&gt;&lt;/li&gt;
&lt;li&gt;&lt;a href=&quot;http://localhost/dotclear/index.php?post/fieldconv&quot;&gt;fieldconv&lt;/a&gt;&lt;/li&gt;
&lt;/ul&gt;
&lt;/div&gt;&lt;/div&gt;</content>
    
    

    
      </entry>
    
  <entry>
    <title>Requirements</title>
    <link href="http://localhost/dotclear/index.php?post/Requirements" rel="alternate" type="text/html"
    title="Requirements" />
    <id>urn:md5:82c57abde809872255772cdbe17d1f3d</id>
    <published>2101-05-08T13:51:00+01:00</published>
    <updated>2013-01-11T04:27:47+01:00</updated>
    <author><name>thierry sousbie</name></author>
        <dc:subject>Install</dc:subject>
            
    <content type="html">    &lt;div&gt;The following libraries/programs are required:&amp;nbsp;&lt;/div&gt;
&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;Cmake 2.6.3+&lt;/li&gt;
&lt;li&gt;GSL&lt;/li&gt;
&lt;/ul&gt;
&lt;br /&gt;
&lt;br /&gt;
&lt;div&gt;The following libraries are optional:&amp;nbsp;&lt;/div&gt;
&lt;br /&gt;
&lt;ul&gt;
&lt;li&gt;CGAL 3.7+ (for delaunay)&lt;/li&gt;
&lt;li&gt;CFitsIO  (for reading FITS)&lt;/li&gt;
&lt;li&gt;SDL/SDL-image (for loading jpg,bmp,...)&lt;/li&gt;
&lt;li&gt;mathGL 1.10 (for visualizing persistence diagrams with pdview)&lt;/li&gt;
&lt;li&gt;Qt4 (for visualizing persistence diagrams with pdview)&lt;/li&gt;
&lt;/ul&gt;
&lt;br /&gt;
&lt;br /&gt;
&lt;div&gt;&lt;ins&gt;NOTE:&lt;/ins&gt; gcc v4.3.3+ is required for multithreading to be enabled.&lt;/div&gt;</content>
    
    

    
      </entry>
  
</feed>