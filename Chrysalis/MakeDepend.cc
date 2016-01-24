// Part of the Makefile system.


#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
//#include <getopt.h>
#include <unistd.h>
#include <algorithm>
#include <functional>
#include <numeric>
#include <iostream>
#include <iterator>
#include <ctype.h>
#include <cstring>
#include <cstdlib>



using namespace std;



inline bool IsSpace(char c) { return isspace(c); }

class recursion_strategy {

 public:
  virtual void GetSubdirs( const string &source_tree_path,
			   set<string> &subdirectories ) = 0;
  virtual ~recursion_strategy() {};
};


class cvs_recursion_strategy : public recursion_strategy {

 public:
  virtual void GetSubdirs( const string &source_tree_path,
			   set<string> &subdirectories );
  virtual ~cvs_recursion_strategy() {};
};


class full_recursion_strategy : public recursion_strategy {

 public:
  virtual void GetSubdirs( const string &source_tree_path,
			   set<string> &subdirectories );
  virtual ~full_recursion_strategy() {};
};



class makefile_builder {

 public:

  makefile_builder( recursion_strategy *rs,
                    bool ignore_archived_directive,
                    bool print_dep_graph,
                    int ubiquitous_cutoff )
    : make_archived( ignore_archived_directive ),
      print_graph( print_dep_graph ),
      ubiq_cutoff( ubiquitous_cutoff ),
      recursion_strat( rs )
  { }

  void ParseDependencies( const string &source_tree_root );

  void GenerateMakefile( const string &makefile );

  void DumpDependencies( const string &filename );

 private:

  class dependency {
    
   public:
    dependency() {}
    dependency( const string &dependent, const string &provider )
      : dependent_(dependent), provider_(provider) { }
    
    const string & Dependent() const { return dependent_; }
    const string & Provider()  const { return provider_; }
    
    bool operator< ( const dependency &other ) const
    { return ( dependent_ < other.dependent_ ||
               dependent_ == other.dependent_ &&
               provider_ < other.provider_ ); }
    
    struct order_by_dependent 
      : public binary_function<dependency,dependency,bool>
    {
      bool operator() ( const dependency &lhs, const dependency &rhs ) const
      { return ( lhs.Dependent() < rhs.Dependent() ); }
    };

    struct has_dependent 
      : public unary_function<dependency,bool>
    {
      has_dependent( const string &dependent ) : dependent_(dependent) {}

      bool operator() ( const dependency &thing ) const
      { return ( thing.Dependent() == dependent_ ); }
     private:
      string dependent_;
    };

   private:
    string dependent_;
    string provider_;
  };

  class subdir_and_target {
   public:
    subdir_and_target() {}
    subdir_and_target( const string &subdir, const string &target ) 
      : subdir_(subdir), target_(target) { }
    
    const string &Subdir() const { return subdir_; }
    const string &Target()   const { return target_; }

    bool operator< ( const subdir_and_target &other ) const
    { return ( subdir_ < other.subdir_ ||
               subdir_ == other.subdir_ &&
               target_ < other.target_ ); }

    bool operator== ( const subdir_and_target &other ) const
    { return ( subdir_ == other.subdir_ &&
               target_ == other.target_ ); }
    
    struct order_by_subdir 
      : public binary_function<subdir_and_target,subdir_and_target,bool>
    {
      bool operator() ( const subdir_and_target &lhs, const subdir_and_target &rhs ) const
      { return ( lhs.Subdir() < rhs.Subdir() ); }
    };

    struct order_by_target 
      : public binary_function<subdir_and_target,subdir_and_target,bool>
    {
      bool operator() ( const subdir_and_target &lhs, const subdir_and_target &rhs ) const
      { return ( lhs.Target() < rhs.Target() ); }
    };
   private:
    string subdir_;
    string target_;
  };

  vector<dependency> compile_time_deps;
  vector<dependency> link_time_deps;
  vector<dependency> run_time_deps;

  set<string> sources_found;
  set<string> headers_found;
  set<string> objects_found;
  set<string> executables_found;

  bool make_archived;
  bool print_graph;
  int ubiq_cutoff;

  recursion_strategy *recursion_strat;

  map<string, string> special_compiles;

  bool IsDirectory( const string &path ) const;
  bool FileIsHeader( const string &filename ) const;
  bool FileIsSource( const string &filename ) const;
  bool FileIsObject( const string &filename ) const;

  void ParseHeaderFile( const string &full_path, const string &partial_path );
  void ParseSourceFile( const string &full_path, const string &partial_path );

  bool CheckLineForInclude( const string &line, string &included_file );
  bool CheckLineForMain( const string &line );
  bool CheckLineForDirective( const string &line, string &directive );

  void RecursiveParseDependencies( const string &source_tree_root,
				   const string &source_tree_subdir );

  void AddHeaderSourceDependencies();

  void GenerateMakefileForObjectFiles( ostream &makefile );

  /// Add to makefile a list of compile commands for .o files whose
  /// .cc files had special compilation directives.
  void GenerateMakefileForSpecialCompiles( ostream &makefile );
  void GenerateMakefileForExecutables( ostream &makefile );

  void GenerateTargetLists( ostream &makefile );

  set<string> nodes;
  set< pair<string,string> > edges;
  
  void PrintGraphNode( ostream& graphStrm, const string& file );
  void PrintGraphEdge( ostream& graphStrm, const string& fromFile, const string& toFile,
                       const bool dashed = false );
};


//
// MAIN
//

int main( int argc, char **argv )
{
  recursion_strategy *rs = new full_recursion_strategy;
  bool ignore_archived_directive = false;
  string target;
  bool print_graph = false;
  int ubiq_cutoff = 0;
  bool print_help = false;

  bool error_in_options = false;
  int option;

  char * options = "Ad:g:u:fh";

  while ( ! error_in_options &&
	  ( ( option = getopt( argc, argv, options ) ) != -1 ) )
    switch ( option ) 
    {
      case 'A' :
        ignore_archived_directive = true;
        break;
      case 'd' :
        target = optarg;
        break;
      case 'f' :
        rs = new full_recursion_strategy;
        break;
      case 'g' :
        target = optarg;
        print_graph = true;
        break;
      case 'u' :
        ubiq_cutoff = atoi( optarg );
        break;
      case 'h' :
        print_help = true;
        break;
      default:
        error_in_options = true;
        print_help = true;
    }

  if ( ! print_graph && ubiq_cutoff > 0 ) {
    cout << "The -u option is only usable in combination with the -g option." << endl << endl;
    print_help = true;
    error_in_options = true;
  }

  if ( print_help )
  {
    cout << "Usage: MakeDepend [-h] [-f] [-d file] [-g file [-u percent]] [-A] dir " << endl
	 << " Parses C and C++ files in dir building a dependencies tree." << endl
	 << "  -h : print this help" << endl
	 << "  -f : recurse over all dirs, not just dirs that are CVS entries" << endl
         << "  -d : dump dependencies of the specified file" << endl
	 << "  -g : print a graph of the dependencies of the specified file in Makefile_auto.<file>.dot" << endl
         << "  -u : the X percent most #included files are excluded from the -g graph" << endl
         << "  -A : ignore 'archived' directives" << endl;
      
    if ( rs )
      delete rs;

    if ( error_in_options )
      return 1;
    else
      return 0;
  }
  

  if ( rs == 0 )
    rs = new cvs_recursion_strategy;

  makefile_builder make_depend( rs, 
                                ignore_archived_directive, 
                                print_graph,
                                ( print_graph ? ubiq_cutoff : 0 ) );

  char * dir = argv[optind];

  make_depend.ParseDependencies( ( dir ? dir : "." ) );


  // If no particular target was specified, use the dependency graph to generate a
  // makefile with the filename "Makefile_auto".  If a particular target was
  // specified, print that target's dependencies to stdout or to a DOT file,
  // depending on the -g option.

  if ( target.empty() )
    make_depend.GenerateMakefile( "Makefile_auto" );
  else
    make_depend.DumpDependencies( target );

  //   cout << "Done!" << endl;

  delete rs;

  return 0;
}



//
// makefile_builder methods
//

// Recursively apply the slave function RecursiveParseDependencies() to
// construct the dependency vectors compile_time_deps, link_time_deps, and
// //run_time_deps from the entire content of the directory subtree rooted at
// source_tree_root.  Each of these three vectors contains ordered pairs (class
// dependency) (A, B) where A depends on B.

void makefile_builder::ParseDependencies( const string &source_tree_root )
{
  this->RecursiveParseDependencies( source_tree_root, "." );

  this->AddHeaderSourceDependencies();

  sort( compile_time_deps.begin(), compile_time_deps.end() );
  sort( link_time_deps.begin(), link_time_deps.end() );
  sort( run_time_deps.begin(), run_time_deps.end() );
}


void makefile_builder::GenerateMakefile( const string &makefile )
{
  ofstream mf( makefile.c_str() );
  
  mf << "###" << "\n";
  mf << "### This file is generated automatically.  DO NOT MODIFY BY HAND!" << "\n";
  mf << "###" << "\n";
  mf << "\n";

  this->GenerateMakefileForObjectFiles( mf );
  this->GenerateMakefileForSpecialCompiles( mf );
  this->GenerateMakefileForExecutables( mf );
  this->GenerateTargetLists( mf );
}
    

void makefile_builder::PrintGraphNode( ostream& graphStrm, const string& file ) {
  string base( file.substr( 0, file.rfind( "." ) ) );
  if ( nodes.count( base ) ) return;
  nodes.insert( base );
  graphStrm << "  \"" << base << "\" [shape=box];" << "\n";
}  


void makefile_builder::PrintGraphEdge( ostream& graphStrm, 
                                       const string& fromFile, const string& toFile,
                                       const bool dashed ) {
  this->PrintGraphNode( graphStrm, fromFile );
  this->PrintGraphNode( graphStrm, toFile );
  string from( fromFile.substr( 0, fromFile.rfind( "." ) ) );
  string   to(   toFile.substr( 0,   toFile.rfind( "." ) ) );
  if ( from == to ) return;
  if ( edges.count( make_pair( from, to ) ) ) return;
  edges.insert( make_pair( from, to ) );
  graphStrm << "  \"" << from << "\" -> \"" << to << "\"";
  if ( dashed ) 
    graphStrm << " [style=dashed]";
  graphStrm << ";" << "\n";
}


void makefile_builder::DumpDependencies( const string &target )
{
  set<dependency> unexpanded_deps;
  set<string> expanded_deps;
  set<string> real_deps, exec_deps, obj_deps;

  ostream* pGraphStrm = 0;
  if ( print_graph ) {
    // Remove slashes and dots from target name so that is safe to use as a
    // filename and put the result in safe_target.
    string safe_target = target;
    string::size_type slashdotpos;
    while ( ( slashdotpos = safe_target.find_first_of( "/." ) ) != string::npos )
      safe_target[slashdotpos] = '_';
    // Open a stream to the file "Makefile_auto.<safe_target>.dot".
    string dotfile = "Makefile_auto.";
    dotfile += safe_target;
    dotfile += ".dot";
    cout << "Saving output to " << dotfile << endl;
    pGraphStrm = new ofstream( dotfile.c_str() );
    *pGraphStrm << "digraph \"" << target << "\" {" << endl;
  }

  // provider_counts[] is, effectively, a symbol table: that is, it's a data
  // structure that maps unique strings to integers, where no string in the
  // table is duplicated.  In the case of this particular symbol table, the
  // strings are the names of code modules (that is, the file names with the
  // extensions such as ".cc" deleted) and the integers are reference counts.
  // If the ordered pair (A, B) is present in the vector link_time_deps - that
  // is, if A depends on B, then the reference count associated with the name B
  // is incremented.

  
  map<string,int> provider_counts;
  for ( vector<dependency>::iterator dep_iter = link_time_deps.begin();
        dep_iter != link_time_deps.end(); ++dep_iter ) {
    string base( dep_iter->Provider().substr( 0, dep_iter->Provider().rfind( "." ) ) );
    provider_counts[base]++;
  }

  for ( vector<dependency>::iterator dep_iter = compile_time_deps.begin();
        dep_iter != compile_time_deps.end(); ++dep_iter ) {
    string base( dep_iter->Provider().substr( 0, dep_iter->Provider().rfind( "." ) ) );
    provider_counts[base]++;
  }

  // The vector of pairs count_per_provider[] simply reverses the map
  // provider_counts[], copying each map element B->n to an ordered pair (n,B),
  // so that the resulting vector can be sorted by number of occurrences, rather
  // than provider name.  (Why not build the mapping as a vector in the first
  // place?  Efficiency: keeping the vector sorted while building up the map
  // would require a lot of copying of data everytime a new provider name was
  // encountered.)

  vector< pair<int,string> > count_per_provider;
  for ( map<string,int>::iterator count_iter = provider_counts.begin();
        count_iter != provider_counts.end(); ++count_iter )
    count_per_provider.push_back( make_pair( count_iter->second, count_iter->first ) );

  // Modules with very high out-degrees - that is, modules on which large 
  // numbers of other modules depend - are defined as "ubiquitous providers" 
  // of dependencies.  More specifically, this class is defined as the top 
  // 'ubiq_cutoff' percent (rounded down) of out-degrees of all dependency
  // providers.

  sort( count_per_provider.begin(), count_per_provider.end() );
  int cutoff_count = count_per_provider.size() * ubiq_cutoff / 100;
  
  // ubiquitous_providers[] is a vector whose elements are the names (strings)
  // of all the providers in this aforementioned class.  It is alphabetically
  // sorted, then dumped to stdout.

  vector<string> ubiquitous_providers;
  for ( int i = 1; i <= cutoff_count; ++i )
    ubiquitous_providers.push_back( count_per_provider[ count_per_provider.size()-i ].second );
  sort( ubiquitous_providers.begin(), ubiquitous_providers.end() );

  if ( ! ubiquitous_providers.empty() ) {
    cout << "Not displaying these \"ubiquitous\" nodes:" << endl;
    copy( ubiquitous_providers.begin(), ubiquitous_providers.end(),
          ostream_iterator<string>( cout, "\n" ) );
  }

  pair<vector<dependency>::iterator,vector<dependency>::iterator> target_range;
  
  // Make the ordered pair target_range the range covering every dependency
  // within link_time_deps[] that specifies a dependency on the given module
  // name 'target'.

  target_range = equal_range( link_time_deps.begin(),
                              link_time_deps.end(),
                              dependency(target,""),
                              dependency::order_by_dependent() );

  unexpanded_deps.insert( target_range.first, target_range.second );
    
  while ( ! unexpanded_deps.empty() )
  {
    dependency a_dep = *(unexpanded_deps.begin());
      
    unexpanded_deps.erase( unexpanded_deps.begin() );
    
    const string &dependent = a_dep.Dependent();
    const string &provider = a_dep.Provider();

    string base( provider.substr( 0, provider.rfind( "." ) ) );
    if ( binary_search( ubiquitous_providers.begin(), ubiquitous_providers.end(), base ) )
      continue;

    // If printing a graph, print only dependencies of a source file on a header
    // file, or dependencies where one or the other file is neither header nor
    // source file.  (This implies that object dependencies will not be printed,
    // and in the unlikely event that there are any header files that depend on
    // source file, such dependencies will not be printed either.)

    if ( print_graph )
      if ( ! this->FileIsSource( dependent ) && ! this->FileIsHeader( dependent ) ||
           ! this->FileIsSource( provider )  && ! this->FileIsHeader( provider ) ||
           this->FileIsHeader( dependent ) && this->FileIsSource( provider ) ) {
        this->PrintGraphEdge( *pGraphStrm, dependent, provider, true );
      }

    // If this dependency provider has not yet been expanded, then insert it
    // into expanded deps, and queue its dependencies for expansion.  If this
    // dependency provider is either a source or a header file, then add it to
    // the set real_deps.  Otherwise, add it to the set obj_deps.
    //
    // Also add any and all run-time dependencies (from run_time_deps) to the
    // set exec_deps, but do not queue these for expansion.

    if ( expanded_deps.count( provider ) == 0 )
    {
      expanded_deps.insert( provider );
      
      pair<vector<dependency>::iterator,vector<dependency>::iterator> more_deps;
	
      more_deps = equal_range( link_time_deps.begin(),
                               link_time_deps.end(),
                               dependency( provider, "" ),
                               dependency::order_by_dependent() );

      for ( ; more_deps.first != more_deps.second; ++more_deps.first )
        unexpanded_deps.insert( *(more_deps.first) );

      if ( this->FileIsSource( provider ) ||
           this->FileIsHeader( provider ) )
        real_deps.insert( provider );
      else
        obj_deps.insert( provider );

      more_deps = equal_range( run_time_deps.begin(),
                               run_time_deps.end(),
                               dependency( provider, "" ),
                               dependency::order_by_dependent() );

      for ( ; more_deps.first != more_deps.second; ++more_deps.first )
        exec_deps.insert( more_deps.first->Provider() );
    }
  }

  if ( this->FileIsObject( target ) )
    obj_deps.insert( target );

  for ( set<string>::iterator obj_dep_iter = obj_deps.begin();
        obj_dep_iter != obj_deps.end(); ++obj_dep_iter ) {
    target_range = equal_range( compile_time_deps.begin(),
                                compile_time_deps.end(),
                                dependency(*obj_dep_iter,""),
                                dependency::order_by_dependent() );
    unexpanded_deps.insert( target_range.first, target_range.second );
  }
  expanded_deps.clear();

  while ( ! unexpanded_deps.empty() )
  {
    dependency a_dep = *(unexpanded_deps.begin());
      
    unexpanded_deps.erase( unexpanded_deps.begin() );
    
    const string &dependent = a_dep.Dependent();
    const string &provider = a_dep.Provider();

    string base( provider.substr( 0, provider.rfind( "." ) ) );
    if ( binary_search( ubiquitous_providers.begin(), ubiquitous_providers.end(), base ) )
      continue;

    if ( print_graph )
      this->PrintGraphEdge( *pGraphStrm, dependent, provider );

    if ( expanded_deps.count( provider ) == 0 )
    {
      expanded_deps.insert( provider );
      
      pair<vector<dependency>::iterator,vector<dependency>::iterator> more_deps;
	
      more_deps = equal_range( compile_time_deps.begin(),
                               compile_time_deps.end(),
                               dependency( provider, "" ),
                               dependency::order_by_dependent() );

      for ( ; more_deps.first != more_deps.second; ++more_deps.first )
        unexpanded_deps.insert( *(more_deps.first) );
    }
  }

  if ( print_graph ) {
    *pGraphStrm << "}" << endl;
    delete pGraphStrm;
  } else {
    copy( real_deps.begin(), real_deps.end(),
          ostream_iterator<string>( cout, "\n" ) );
    copy( exec_deps.begin(), exec_deps.end(),
          ostream_iterator<string>( cout, "\n" ) );
  }
}

  
// This function is effectively a mapcar of the subordinate functions
// Parse{Source,Header}File over the content of source_tree_subdir.  These two
// subordinate functions add dependencies, represented as ordered pairs (A, B)
// where A depends on B, to the three vectors of class dependency
// compile_time_deps, link_time_deps, and run_time_deps.

void makefile_builder::RecursiveParseDependencies( const string &source_tree_root, 
						   const string &source_tree_subdir )
{
  string subdir_path = source_tree_root + "/" + source_tree_subdir;

  const char* dirname = subdir_path.c_str();
    
  struct stat dir_stat;
  if ( ( stat( dirname, &dir_stat ) == 0 ) &&
       ( S_ISDIR( dir_stat.st_mode ) ) ) 
  {
    if ( DIR* dir_ptr = opendir( dirname ) ) 
    {
      while (struct dirent* dir_entry = readdir(dir_ptr)) 
      {
	string entryname( dir_entry->d_name );
	  
	if ( entryname[0] == '.' )
	  continue;
	
	if ( "MakeDepend.cc" == entryname )
	  continue;

	string partial_path = source_tree_subdir + "/" + entryname;

	while ( partial_path.substr( 0, 2 ) == "./" )
	  partial_path.erase( 0, 2 );

	string full_path = source_tree_root + "/" + partial_path;

	if ( this->FileIsHeader( entryname ) )
	{
	  // 	  cout << "Parsing header file " << partial_path << "." << endl;
	  this->ParseHeaderFile( full_path, partial_path );
	}
	else if ( this->FileIsSource( entryname ) )
	{
	  // 	  cout << "Parsing source file " << partial_path << "." << endl;
	  this->ParseSourceFile( full_path, partial_path );
	}
      }
    }
  }

  set<string> subdirs;
  recursion_strat->GetSubdirs( subdir_path, subdirs );

  //   if ( ! subdirs.empty() ) {
  //     cout << "Recursing over:" << endl;
  //     copy( subdirs.begin(), subdirs.end(),
  //           ostream_iterator<string>( cout, "\n" ) );
  //   }

  set<string>::iterator subdir_iter;
  for ( subdir_iter = subdirs.begin(); subdir_iter != subdirs.end(); ++subdir_iter )
    this->RecursiveParseDependencies( source_tree_root, 
				      source_tree_subdir + "/" + *subdir_iter );
}



void cvs_recursion_strategy::GetSubdirs( const string &source_tree_dir,
					 set<string> &subdirectories )
{
  vector< pair<string,string> > file_and_prefix;
  file_and_prefix.push_back( make_pair( source_tree_dir + "/CVS/Entries", "D/" ) );
  file_and_prefix.push_back( make_pair( source_tree_dir + "/CVS/Entries.Log", "A D/" ) );

  for ( unsigned int i = 0; i < file_and_prefix.size(); ++i ) {
    string entries_file = file_and_prefix[i].first;
    string dir_entry_prefix = file_and_prefix[i].second;
    struct stat file_stat;
    if ( ( stat( entries_file.c_str(), &file_stat ) == 0 ) &&
         ( S_ISREG( file_stat.st_mode ) ) ) 
    {
      ifstream entries( entries_file.c_str() );

      const int bufsize = 8092;
      char buffer[bufsize];
      string line;
    
      while ( entries.getline( buffer, bufsize ) )
      {
        line = buffer;
        if ( 0 == line.compare( 0, dir_entry_prefix.size(), dir_entry_prefix ) && 
             line.find( '/' ) != line.rfind( '/' ) )
        {
          line = line.erase( 0, dir_entry_prefix.size() );
          subdirectories.insert( line.substr( 0, line.find( '/') ) );
        }
      }
    }
  }
}
  


void full_recursion_strategy::GetSubdirs( const string &source_tree_dir,
					  set<string> &subdirectories )
{
  const char* dirname = source_tree_dir.c_str();
    
  if ( DIR* dir_ptr = opendir( dirname ) ) 
  {
    string dir_entry_path;
    struct stat dir_stat;

    while (struct dirent* dir_entry = readdir(dir_ptr)) 
    {
      dir_entry_path = source_tree_dir + "/" + dir_entry->d_name;
      
      if ( ( stat( dir_entry_path.c_str(), &dir_stat ) == 0 ) &&
	   ( S_ISDIR( dir_stat.st_mode ) ) )
      {
	if ( dir_entry->d_name[0] != '.' ) 
	  subdirectories.insert( dir_entry->d_name );
      }
    }
  }
}



void makefile_builder::AddHeaderSourceDependencies()
{
  for ( set<string>::iterator source_iter = sources_found.begin();
	source_iter != sources_found.end(); ++source_iter )
  {
    if ( ! this->FileIsSource( *source_iter ) )
      continue;

    const string &source_file = *source_iter;

    string source_base( source_file.substr( 0, source_file.rfind( "." ) ) );
    string corresp_header1( source_base + ".h" );
    string corresp_header2( source_base + ".hh" );

    if ( headers_found.count( corresp_header1 ) )
      link_time_deps.push_back( dependency( corresp_header1, source_file ) );
		     
    if ( headers_found.count( corresp_header2 ) )
      link_time_deps.push_back( dependency( corresp_header2, source_file ) );
  }
}	     
  


void makefile_builder::GenerateMakefileForObjectFiles( ostream &mf )
{
  pair<vector<dependency>::iterator,vector<dependency>::iterator> compile_target;
  
  compile_target.second = compile_time_deps.begin();
  
  while ( compile_target.second != compile_time_deps.end() )
  {
    compile_target = equal_range( compile_time_deps.begin(),
				  compile_time_deps.end(),
				  *(compile_target.second),
				  dependency::order_by_dependent() );

    if ( ! this->FileIsObject( compile_target.first->Dependent() ) )
      continue;

    //     cout << "Building makefile for compile target: " 
    //  	 << compile_target.first->Dependent()
    //  	 << endl;

    set<dependency> unexpanded_deps;
    set<string> expanded_deps;
    set<string> real_deps;

    string target = compile_target.first->Dependent();

    unexpanded_deps.insert( compile_target.first, compile_target.second );
    
    while ( ! unexpanded_deps.empty() )
    {
      dependency a_dep = *(unexpanded_deps.begin());
      
      //       cout << "Looking at " << a_dep.Dependent() << " -> " << a_dep.Provider() << endl;

      unexpanded_deps.erase( unexpanded_deps.begin() );

      if ( expanded_deps.count( a_dep.Provider() ) == 0 )
      {
	expanded_deps.insert( a_dep.Provider() );

	pair<vector<dependency>::iterator,vector<dependency>::iterator> more_deps;
	
	more_deps = equal_range( compile_time_deps.begin(),
				 compile_time_deps.end(),
				 dependency( a_dep.Provider(), "" ),
				 dependency::order_by_dependent() );

	//  	cout << "Found " << distance( more_deps.first, more_deps.second ) 
	//  	     << " more deps, " << flush;

	// 	int before = unexpanded_deps.size();
	for ( ; more_deps.first != more_deps.second; ++more_deps.first )
	  if ( expanded_deps.count( more_deps.first->Provider() ) == 0 )
	    unexpanded_deps.insert( *(more_deps.first) );
	// 	int after = unexpanded_deps.size();

	//  	cout << after - before << " of which are new." << endl;
	

	real_deps.insert( "$(SRC)/" + a_dep.Provider() );
	// Is file there???
	//cout << "Testing: " << a_dep.Provider().c_str() << endl;
	//FILE * pTest = fopen(a_dep.Provider().c_str(), "r");
	//if (pTest != NULL)
	//  real_deps.insert( "$(SRC)/" + a_dep.Provider() );
	//else 
	//  real_deps.insert( "$(SRC)/Spines/" + a_dep.Provider() );
	
	//if (pTest != NULL)
	//  fclose(pTest);
      }
    }

    mf << target << ": ";
    copy( real_deps.begin(), real_deps.end(),
	  ostream_iterator<string>( mf, " " ) );
    mf << "\n";
  }
}

void makefile_builder::GenerateMakefileForSpecialCompiles( ostream &mf )
{
  map<string,string>::iterator iter;
  string ccfile;
  for (iter = special_compiles.begin();iter != special_compiles.end();++iter){
    //ccfile = iter->first.substr(0, iter->first.rfind('.')) + ".cc";
    mf << iter->first << ":\n"
       << "\t@ mkdir -p $(OBJ)/${@D}\n"
       << "\t$(CPLUSPLUS) $(CPPC) " << iter->second 
       << " -c $(SRC)/$*.cc -o $(OBJ)/$*.o\n";
  }
}




void makefile_builder::GenerateMakefileForExecutables( ostream &mf )
{
  pair<vector<dependency>::iterator,vector<dependency>::iterator> link_target;
  
  link_target.second = link_time_deps.begin();
  
  while ( link_target.second != link_time_deps.end() )
  {
    link_target = equal_range( link_time_deps.begin(),
			       link_time_deps.end(),
			       *(link_target.second),
			       dependency::order_by_dependent() );

    //     cout << "On " << link_target.first->Dependent() << ".  " << flush;
    //     cout << distance( link_target.second, link_time_deps.end() ) << " to go." << endl;

    if ( this->FileIsHeader( link_target.first->Dependent() ) ||
	 this->FileIsSource( link_target.first->Dependent() ) ||
	 this->FileIsObject( link_target.first->Dependent() ) )
      continue;

    //     cout << "Building makefile for link target: " 
    //   	 << link_target.first->Dependent()
    //   	 << endl;

    set<dependency> unexpanded_deps;
    set<string> expanded_deps;
    set<string> obj_deps;
    set<string> exec_deps;

    bool needs_xerces_lib = false;

    unexpanded_deps.insert( link_target.first, link_target.second );
    
    while ( ! unexpanded_deps.empty() )
    {
      dependency a_dep = *(unexpanded_deps.begin());
      
      //       cout << "Looking at " << a_dep.Dependent() << " -> " << a_dep.Provider() << endl;

      if ( a_dep.Provider().find( "xerces" ) != string::npos )
	needs_xerces_lib = true;

      pair<vector<dependency>::iterator,vector<dependency>::iterator> rt_deps;
      rt_deps = equal_range( run_time_deps.begin(),
                             run_time_deps.end(),
                             dependency( a_dep.Provider(), "" ),
                             dependency::order_by_dependent() );
      for ( ; rt_deps.first != rt_deps.second; ++rt_deps.first )
        exec_deps.insert( rt_deps.first->Provider() );

      unexpanded_deps.erase( unexpanded_deps.begin() );

      if ( expanded_deps.count( a_dep.Provider() ) == 0 )
      {
	expanded_deps.insert( a_dep.Provider() );

	pair<vector<dependency>::iterator,vector<dependency>::iterator> more_deps;
	more_deps = equal_range( link_time_deps.begin(),
				 link_time_deps.end(),
				 dependency( a_dep.Provider(), "" ),
				 dependency::order_by_dependent() );

	//  	cout << "Found " << distance( more_deps.first, more_deps.second ) 
	//  	     << " more deps, " << flush;

	// 	int before = unexpanded_deps.size();
	for ( ; more_deps.first != more_deps.second; ++more_deps.first )
	  if ( expanded_deps.count( more_deps.first->Provider() ) == 0 )
	    unexpanded_deps.insert( *(more_deps.first) );
	// 	int after = unexpanded_deps.size();
	
	//  	cout << after - before << " of which are new." << endl;

	if ( this->FileIsObject( a_dep.Provider() ) )
	{
	  obj_deps.insert( a_dep.Provider() );
	}
      }
    }


    // Executable dependencies
    mf << link_target.first->Dependent() << ": ";
    for ( set<string>::iterator dep_iter = exec_deps.begin(); 
	  dep_iter != exec_deps.end(); ++dep_iter )
      mf << " " << *dep_iter;
    mf << " checkLock";
    mf << "\n";

    // Object dependencies
    mf << link_target.first->Dependent() << ": ";
    for ( set<string>::iterator dep_iter = obj_deps.begin(); 
	  dep_iter != obj_deps.end(); ++dep_iter )
      mf << " " << *dep_iter;
    mf << "\n";

    // Generate instructions to build main program executable.  The second for
    // loop does nothing more than select the main program.

    string targetname = link_target.first->Dependent( );
    string libname = "_" + targetname + "_temp";
    mf << "\t" << "/bin/rm -f lib" << libname << ".a\n";
    mf << "\t" << "$(AR) $(ARFLAGS) lib" << libname << ".a";
    for ( set<string>::iterator dep_iter = obj_deps.begin(); 
          dep_iter != obj_deps.end(); ++dep_iter )
    {    
      mf << " $(OBJ)/" << *dep_iter;     
    }
    mf << "\n\t$(BIN)/checkLock $(BIN)/$@";
    mf << "\n\t$(CPLUSPLUS) $(LINK_OPTIONS) -o $(BIN)/$@ ";
    for ( set<string>::iterator dep_iter = obj_deps.begin();
          dep_iter != obj_deps.end(); ++dep_iter )
    {
      string tofind = "/" + targetname + ".o";
      string context = "/" + *dep_iter;
      string::size_type found = context.rfind(tofind);
      if ( found == string::npos ) continue;
      if ( found + tofind.size( ) == context.size( ) )
      {    
        mf << " $(OBJ)/" << *dep_iter;
      }
    }
    mf << " -L. $(LINK_LIBS) -l" << libname;
    if (needs_xerces_lib) mf << " $(XERCES_LIB)";
    mf << "\n\t" << "/bin/rm lib" << libname << ".a\n";
  }
}    


void makefile_builder::GenerateTargetLists( ostream &mf )
{
  mf << "\n";
  mf << "# Target lists\n";

  mf << "\n";
  mf << "OBJECTS = \\\n";

  set<string>::iterator object_iter = objects_found.begin();
  for ( ; object_iter != objects_found.end(); ++object_iter )
    mf << "  " << *object_iter << "\\\n";
  
  mf << "\n";
  mf << "nolink: $(OBJECTS)\n";

  mf << "\n";
  mf << "EXECUTABLES = \\\n";

  vector<subdir_and_target> targets_by_subdir;

  vector<dependency>::iterator link_target = link_time_deps.begin();
  for ( ; link_target != link_time_deps.end(); ++link_target )
  {
    if ( this->FileIsHeader( link_target->Dependent() ) ||
	 this->FileIsSource( link_target->Dependent() ) ||
	 this->FileIsObject( link_target->Dependent() ) )
      continue;

    const string &provider = link_target->Provider();

    // grab the directory part of the provider's path
    string::size_type slash_pos = provider.rfind('/');
    if ( slash_pos == string::npos )
      slash_pos = 0;
    string subdir = provider.substr( 0, slash_pos );
    
    // change all slashes to dots while simultaneously making each
    // subdir a dependency of its parent so that making 'foo' will not
    // only make the executables in 'foo' but also those in 'foo/bar'.
    string parentdir, currdir;
    while ( ( slash_pos = subdir.find('/') ) != string::npos )
    {
      currdir = subdir.substr( 0, slash_pos );
      if ( ! parentdir.empty() )
        targets_by_subdir.push_back( subdir_and_target( parentdir, currdir ) );
      parentdir = currdir;

      subdir[slash_pos] = '.';
    }
    if ( ! parentdir.empty() )
      targets_by_subdir.push_back( subdir_and_target( parentdir, subdir ) );
    
    if ( ! subdir.empty() )
      targets_by_subdir.push_back( subdir_and_target( subdir, link_target->Dependent() ) );

    mf << "  " << link_target->Dependent() << "\\\n";
  }
  
  mf << "\n";
  mf << "everything: $(EXECUTABLES)\n\n";

  sort( targets_by_subdir.begin(), targets_by_subdir.end() );
  targets_by_subdir.erase( unique( targets_by_subdir.begin(), 
                                   targets_by_subdir.end() ),
                           targets_by_subdir.end() );

  pair<vector<subdir_and_target>::iterator,vector<subdir_and_target>::iterator> range;
  range.second = targets_by_subdir.begin();
  while ( range.second != targets_by_subdir.end() )
  {
    range = equal_range( targets_by_subdir.begin(),
                         targets_by_subdir.end(),
                         *(range.second),
                         subdir_and_target::order_by_subdir() );

    const string &subdir = range.first->Subdir();
    mf << ".PHONY: " << subdir << "\n";
    mf << subdir << ": \\\n";
    
    for ( ; range.first != range.second; ++range.first )
      mf << "  " << range.first->Target() << "\\\n";
    mf << "\n";
  }
}


void makefile_builder::ParseHeaderFile( const string &full_path, const string &partial_path )
{
  const string &header_file = partial_path;

  headers_found.insert( header_file );

  ifstream header( full_path.c_str() );
  
  if ( ! header )
  {
    cerr << "Unable to open '" << full_path << "' for reading." << endl;
    cerr << "Exiting." << endl;
    exit( -1 );
  }

  const int bufsize = 8092;
  char buffer[bufsize];
  string line;
  string included_file;
  string directive;

  while ( header.getline( buffer, bufsize ) )
  {
    line = buffer;
    
    if ( this->CheckLineForInclude( line, included_file ) )
    {
      string prefix;
      if (strstr(header_file.c_str(), "Spines") != NULL)
	prefix = "Spines/";

      dependency new_dep( header_file, prefix + included_file );
     
      //if (strstr(header_file.c_str(), "System") != NULL)
      //cerr << "Included: " << included_file <<" from " << header_file <<  endl;
      compile_time_deps.push_back( new_dep );
      link_time_deps.push_back( new_dep );
    }
    else if ( this->CheckLineForDirective( line, directive ) ) {
      cerr << "Found a MakeDepend directive in the file" << endl;
      cerr << header_file << endl;
      cerr << "Directives are not allowed in header files." << endl;
      exit( -1 );
    }
  }
}
 

void makefile_builder::ParseSourceFile( const string &full_path, const string &partial_path )
{
  const string &source_file = partial_path;

  sources_found.insert( source_file );
  
  ifstream source( full_path.c_str() );
  
  if ( ! source )
  {
    cerr << "Unable to open '" << full_path << "' for reading." << endl;
    cerr << "Exiting." << endl;
    exit( -1 );
  }

  string sourcefile_base( source_file.substr( 0, source_file.rfind('.') ) );
  string object_file = sourcefile_base + ".o";

  string filename( source_file.substr( source_file.rfind('/') + 1 ) );
  string filename_base( filename.substr( 0, filename.rfind('.') ) );
  string executable = filename_base;

  const int bufsize = 8092;
  char buffer[bufsize];
  string line;
  string included_file;
  bool is_main = false;
  bool found_include = false;
  bool bSkip = false;

  if (strstr(source_file.c_str(), "Spines") != NULL)
    bSkip = true;


  string directive;
  string compflag;
  vector<string> explicit_dependencies;

  while ( source.getline( buffer, bufsize ) ) {
    line = buffer;

    if ( this->CheckLineForDirective( line, directive ) ) {
      if ( directive.substr(0,12) == " dependency " ) {
        explicit_dependencies.push_back( directive.substr(12) );
      } 
      else if ( directive == " archived" ) {    
	if ( ! make_archived ) return;   
      }
      else if (directive.substr(0,4) == " CC_" || directive.substr(0, 7) == " cflags") {
	if (directive.substr(0,4) == " CC_")
	  compflag = directive.substr(4);
	else
	  compflag = directive.substr(7);
	  
	//remove all whitespace
	compflag.erase(remove_if(compflag.begin(),compflag.end(),
				  IsSpace), 
			compflag.end());
	compflag = "$(" + compflag + ") ";
	//either create the flag or just append to it; 
        special_compiles[object_file] += compflag;
      }
      else {
	if (directive.substr(0,8) != " library" && line != "// MakeDepend: shared SortKmers") {
	  cerr << "Warning: Unknown MakeDepend directive: ";
	  cerr << line;
	  cerr << " in " << full_path << ". ";
	  cerr << "Continuing...\n";
	}
	//exit(1);
      }
    }

    // Check for main program.
    else if ( !bSkip && this->CheckLineForMain( line ) ) {
      is_main = true;
    }
    else if ( this->CheckLineForInclude( line, included_file ) ) {
      found_include = true;
      string prefix;
      if (strstr(source_file.c_str(), "Spines") != NULL) {
	prefix = "Spines/";
      }
      //cerr << "Using prefix for file " << included_file << endl;
	//}

      compile_time_deps.push_back( dependency( source_file, prefix + included_file ));
      link_time_deps.push_back( dependency( source_file, prefix + included_file ) );
    }
  }

  objects_found.insert( object_file );
  compile_time_deps.push_back( dependency( object_file, source_file ) );

  if ( is_main ) {
    if ( executables_found.find( executable ) != executables_found.end() ) {
      // Conflicting executables.

      vector<dependency>::iterator link_target;
  
      // Find other object file.
      link_target = find_if( link_time_deps.begin(),
			     link_time_deps.end(),
			     dependency::has_dependent( executable ) );
      
      if ( link_target == link_time_deps.end() ) {
        cerr << "Unknown error caused by executable target collision." << endl;
        cerr << "There are two source files that will produce identically named" << endl;
        cerr << "executables.  Until this is resolved, compilation cannot continue." << endl;
        exit(1);
      }

      string other_object_file = link_target->Provider();

      link_target = find_if( link_time_deps.begin(), 
                             link_time_deps.end(),
                             dependency::has_dependent( other_object_file ) );

      if ( link_target == link_time_deps.end() )
	{
	  cerr << "Unknown error caused by executable target collision." << endl;
	  cerr << "There are two source files that will produce identically named" << endl;
	  cerr << "executables.  Until this is resolved, compilation cannot continue." << endl;
	  exit(1);
	}

      string other_source_file = link_target->Provider();

      cerr << "There is a conflict between two main programs." << endl;
      cerr << "The executable '" << executable << "'"
           << " could be built from either of these files: " << endl;
      cerr << "  " << source_file << endl;
      cerr << "  " << other_source_file << endl;
      cerr << "Until one of these source files is removed or renamed,"
           << " compilation cannot continue." << endl;
      exit(1);
    }

    executables_found.insert( executable );
    
    // If this is a main program, when linking we'll need to be able
    // to get from the executable to the object to the source (and
    // thence to the other headers/sources/objects).
    link_time_deps.push_back( dependency( executable, object_file ) );
    link_time_deps.push_back( dependency( object_file, source_file ) );
  }
  else
  {
    // If this is not a main program, when linking we'll need to be
    // able to get from the source to object file, so executables can
    // find the object file.
    link_time_deps.push_back( dependency( source_file, object_file ) );
  }

  // Add in the explicit dependencies.
  for ( vector<string>::iterator exp_iter = explicit_dependencies.begin();
        exp_iter != explicit_dependencies.end(); ++exp_iter )
    run_time_deps.push_back( dependency( object_file, *exp_iter ) );
}
 

bool makefile_builder::CheckLineForDirective( const string &line, string& directive )
{
  directive.clear();

  const unsigned int DIRECTIVE_START=14;
  if ( line.size( ) >= DIRECTIVE_START && line[0] == '/'
       && strncmp( line.c_str(), "// MakeDepend:", DIRECTIVE_START) == 0) 
    directive = line.substr(DIRECTIVE_START);

  return ( ! directive.empty() );
}


bool makefile_builder::CheckLineForInclude( const string &line, string &included_file )
{
  if ( line.size() < 8 )
    return false;

  if ( line[0] != '#' ||
       strncmp( line.c_str(), "#include", 8 ) != 0 ) 
    return false;

  string::size_type first_quote = line.find( '"', 8 );
  
  if ( first_quote == string::npos ) // if not found
    return false;
  
  string::size_type first_letter = first_quote + 1;
  string::size_type next_quote = line.find( '"', first_letter );
  
  if ( next_quote == string::npos ) // if not found
    return false;
  
  included_file = line.substr( first_letter, next_quote-first_letter );
  return true;
}


bool makefile_builder::CheckLineForMain( const string &line )
{
  return ( strncmp( line.c_str(), "int main(", 9 ) == 0 );
}


bool makefile_builder::IsDirectory( const string &path ) const
{
  struct stat path_stat;
  if ( stat( path.c_str(), &path_stat ) != 0 )
  {
    cerr << "stat() failed on path '" << path << "'." << endl;
    cerr << "Exiting." << endl;
    exit( -1 );
  }

  return ( S_ISDIR( path_stat.st_mode ) );
}


int CompareSubstring( const string &the_string, 
                      unsigned int start, 
                      unsigned int length, 
                      const char *query )
{
#if __GNUC__ > 2
  return the_string.compare( start, length, query );
#else
  return the_string.compare( query, start, length );
#endif
}


bool makefile_builder::FileIsHeader( const string &filename ) const
{
  return ( filename.size() > 2 &&
	   ( CompareSubstring( filename, filename.size() - 2, 2, ".h" ) == 0 ) ||
	   filename.size() > 3 &&
	   ( CompareSubstring( filename, filename.size() - 3, 3, ".hh" ) == 0 ) );
}


bool makefile_builder::FileIsSource( const string &filename ) const
{
  return ( filename.size() > 3 &&
	   CompareSubstring( filename, filename.size() - 3, 3, ".cc" ) == 0 );
}


bool makefile_builder::FileIsObject( const string &filename ) const
{
  return ( filename.size() > 2 &&
	   CompareSubstring( filename, filename.size() - 2, 2, ".o" ) == 0 );
}


