#!/usr/bin/perl

#======================================== DEPENDENCIES =========================================================================================
use lib "$ENV{GUS_HOME}/lib/perl";
use strict;
use Data::Dumper;
use Getopt::Long;
use GUS::ObjRelP::DbiDatabase;
use GUS::Supported::GusConfig;
use DBI;
use DBD::Oracle;

#===================================== ARGUMENTS AND SETUP ====================================================================================
my ($testFile, $gusConfigFile, $sequenceId, $regionFile);

&GetOptions("test_file=s"=> \$testFile,
            "sequence_id=s"=> \$sequenceId,
           );

$gusConfigFile = $ENV{GUS_HOME} . "/config/gus.config";
my $gusConfig = GUS::Supported::GusConfig->new($gusConfigFile);

open(TEST, ">$testFile") or die "Cannot create a test file: $!";

my $db = GUS::ObjRelP::DbiDatabase->new($gusConfig->getDbiDsn(),
                                         $gusConfig->getDatabaseLogin(),
                                         $gusConfig->getDatabasePassword(),
                                         0, 0, 1,
                                         $gusConfig->getCoreSchemaName()
    );

my $dbh = $db->getQueryHandle();

#======================================= GENERATING LOCATIONSHIFTS ARRAY  =========================================================================

my $INDEL_QUERY = "select i.location, i.shift
                      from apidb.indel i
                      where sample_name = '$sequenceId'";

my $INDEL_QUERY_SH = $dbh->prepare($INDEL_QUERY);

$INDEL_QUERY_SH->execute();

my @locationshifts = ();
my $counter = 0;
my $currentShift = 0;
while (my ($location, $shift) = $INDEL_QUERY_SH->fetchrow_array()) {
    push ( @{$locationshifts[$counter]}, ($location, $shift + $currentShift));
    $counter++;
    $currentShift = $shift + $currentShift;
}
#print TEST Dumper(\@locationshifts), "\n";
#close TEST;

#====================================== GENERATING COORDINATES ARRAY ============================================================================

my $COORDINATE_QUERY = "select listagg(taf.source_id, ',') WITHIN GROUP (ORDER BY taf.aa_feature_id) as transcripts, s.source_id, tf.parent_id as gene_na_feature_id, el.start_min as exon_start, el.end_max as exon_end, decode(el.is_reversed, 1, afe.coding_end, afe.coding_start) as cds_start, decode(el.is_reversed, 1, afe.coding_start, afe.coding_end) as cds_end, el.is_reversed from dots.transcript tf , dots.translatedaafeature taf , dots.aafeatureexon afe , dots.exonfeature ef , dots.nalocation el, dots.nasequence s where tf.na_feature_id = taf.na_feature_id and taf.aa_feature_id = afe.aa_feature_id and afe.exon_feature_id = ef.na_feature_id and ef.na_feature_id = el.na_feature_id  and ef.na_sequence_id = s.na_sequence_id and tf.external_database_release_id = 5981 and s.source_id = 'OU755535' group by s.source_id, tf.parent_id, el.start_min, el.end_max, afe.coding_start, afe.coding_end, el.is_reversed order by s.source_id, el.start_min";
my $COORDINATE_QUERY_SH = $dbh->prepare($COORDINATE_QUERY);
$COORDINATE_QUERY_SH->execute();
my @coordinates = ();
$counter=0;
while (my ($transcript, $source_id, $gene_na_feature_id, $exon_start, $exon_end, $cds_start, $cds_end, $rev) = $COORDINATE_QUERY_SH->fetchrow_array()) {
  push (@{$coordinates[$counter]}, ($cds_start, $cds_end, $rev));
  $counter++;
}
#print TEST Dumper(\@coordinates), "\n";
#close TEST;

#========================== Generating coordinatesNonCodingLen Object ===========================================================================

my @coordinatesNonCodingLen = ();
my $noncodingLen=0;
my $endingPosition=0;
my $coordinatesLen = scalar @coordinates;
for (my $i=0;$i<$coordinatesLen;$i++){
    $noncodingLen = $coordinates[$i][0] - $endingPosition - 1 + $noncodingLen;
    $endingPosition = $coordinates[$i][1];
    push (@{$coordinatesNonCodingLen[$i]}, ($coordinates[$i][0], $coordinates[$i][1], $noncodingLen));
}
#print TEST Dumper(\@coordinatesNonCodingLen), "\n";
#close TEST;

#======================================= GENERATING TRANSCRIPT LOCATION ARRAY  =========================================================================                                                    
my @transcriptLocations = ();
my $coordinatesNonCodingFrame=0;
my $shiftFrame=0;
my $coordinatesNonCodingFrameLen = scalar @coordinatesNonCodingLen;
my $coordinatesNonCodingFrameLimit = $coordinatesNonCodingFrameLen - 1;
my $locationShiftsLen = scalar @locationshifts;
my $snpLocation;
my $transcriptLocation;
my $lastCodingShiftFrame=0;
for (my $shiftFrame;$shiftFrame<$locationShiftsLen;$shiftFrame++){
    $snpLocation = $locationshifts[$shiftFrame][0];
    ($transcriptLocation, $coordinatesNonCodingFrame,$lastCodingShiftFrame) = &getTranscriptLocation($snpLocation, $coordinatesNonCodingFrame, $coordinatesNonCodingFrameLimit, $shiftFrame, $lastCodingShiftFrame);
    push (@{$transcriptLocations[$shiftFrame]}, ($snpLocation, $transcriptLocation));
}
print TEST Dumper(\@transcriptLocations), "\n";
close TEST; 

#=================================== SUBROUTINES ==========================================================================================

sub getTranscriptLocation {
    my ($snpLocation, $coordinatesNonCodingFrame, $coordinatesNonCodingFrameLimit, $shiftFrame, $lastCodingShiftFrame) = @_;
    my $transcriptLocation;
    if ($coordinatesNonCodingFrame == $coordinatesNonCodingFrameLimit) {
        if ($snpLocation >= $coordinatesNonCodingLen[$coordinatesNonCodingFrame][0] && $snpLocation <= $coordinatesNonCodingLen[$coordinatesNonCodingFrame][1]) {
            ($transcriptLocation, $lastCodingShiftFrame) = &calcTranscriptLocation($snpLocation, $shiftFrame, $coordinatesNonCodingFrame, $lastCodingShiftFrame);
        }
        else {
	    $transcriptLocation = "NC";
        }
    }
    elsif ($snpLocation >= $coordinatesNonCodingLen[$coordinatesNonCodingFrame][0] && $snpLocation <= $coordinatesNonCodingLen[$coordinatesNonCodingFrame][1]) {
	($transcriptLocation, $lastCodingShiftFrame) = &calcTranscriptLocation($snpLocation, $shiftFrame, $coordinatesNonCodingFrame, $lastCodingShiftFrame);
    }
    elsif ($snpLocation > $coordinatesNonCodingLen[$coordinatesNonCodingFrame][1] || $coordinatesNonCodingFrame == $coordinatesNonCodingFrameLimit) {
        until ($snpLocation <= $coordinatesNonCodingLen[$coordinatesNonCodingFrame][1] || $coordinatesNonCodingFrame == $coordinatesNonCodingFrameLimit) {
            $coordinatesNonCodingFrame++;
        }
        if ($coordinatesNonCodingFrame == $coordinatesNonCodingFrameLimit) {
            if ($snpLocation >= $coordinatesNonCodingLen[$coordinatesNonCodingFrame][0] && $snpLocation <= $coordinatesNonCodingLen[$coordinatesNonCodingFrame][1]) {
		($transcriptLocation, $lastCodingShiftFrame) = &calcTranscriptLocation($snpLocation, $shiftFrame, $coordinatesNonCodingFrame, $lastCodingShiftFrame);
            }
            else {
                $transcriptLocation = "NC";
            }
        }
        elsif ($snpLocation >= $coordinatesNonCodingLen[$coordinatesNonCodingFrame][0]) {
	    ($transcriptLocation, $lastCodingShiftFrame) = &calcTranscriptLocation($snpLocation, $shiftFrame, $coordinatesNonCodingFrame, $lastCodingShiftFrame);
        }
        else {
            $transcriptLocation = "NC";
        }
    }
    else {
        $transcriptLocation = "NC";
    }
    return ($transcriptLocation, $coordinatesNonCodingFrame, $lastCodingShiftFrame);
}


sub calcTranscriptLocation {
    my ($snpLocation, $shiftFrame, $coordinatesNonCodingFrame, $lastCodingShiftFrame) = @_;
    my $transcriptLocation;
    if ($shiftFrame == 0) {
        $transcriptLocation = $snpLocation - $coordinatesNonCodingLen[$coordinatesNonCodingFrame][2];
	$lastCodingShiftFrame = $shiftFrame;
    }
    else {
        $transcriptLocation = $snpLocation - $coordinatesNonCodingLen[$coordinatesNonCodingFrame][2] + $locationshifts[$lastCodingShiftFrame][1];
	$lastCodingShiftFrame = $shiftFrame;
    }
    return ($transcriptLocation, $lastCodingShiftFrame);
}

