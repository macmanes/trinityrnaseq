#!/usr/bin/env perl

use strict;
use warnings;

use DBI;

my $obo_tab_file = $ARGV[0] or die "usage: $0 obo.tab\n\n";

main: {

    if (-s "GO.sqlite") {
        die "Error, remove existing GO.sqlite database before regenerating it here.\n";
    }
    

    my $dbh = DBI->connect( "dbi:SQLite:GO.sqlite" ) || die "Cannot connect: $DBI::errstr";

    $dbh->do("PRAGMA synchronous=OFF");

    my $create_GO_table_ddl = qq {
        CREATE TABLE go (
                         id varchar(20),
                         name TEXT,
                         namespace varchar(30),
                         def TEXT
        )
        };

    $dbh->do($create_GO_table_ddl);

    $dbh->do("CREATE UNIQUE INDEX id_idx ON go(id)");
    
    my $insert_entry_dml = qq {
        INSERT INTO go
        VALUES (?,?,?,?)
    };

    my $insert_entry_dsh = $dbh->prepare($insert_entry_dml);


    my $counter = 0;
    open (my $fh, $obo_tab_file) or die $!;
    while (<$fh>) {
        chomp;
        my ($id, $name, $namespace, $def) = split(/\t/);
        $insert_entry_dsh->execute($id, $name, $namespace, $def);
        
        $counter++;
        print STDERR "\r[$counter]   ";
    }
    close $fh;

    $insert_entry_dsh->finish();
    
    print STDERR "\n\ndone.\n\n";
    

    $dbh->disconnect();

    exit(0);
}
