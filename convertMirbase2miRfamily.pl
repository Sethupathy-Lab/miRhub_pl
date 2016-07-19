#!/usr/bin/perl -w

$miRlist=shift;
open (LIST,$miRlist) or die;
while(<LIST>) {
   chomp();
   $A0{$_} =1;
}
close(LIST);

$miRcoords='hsa_tableL.bed';
$FamFile='miR_Family_Info.txt';
open (AFFY,$FamFile) or die;
$head=<AFFY>;
while(<AFFY>) {
   chomp();

   my($fam,$seed,$sp,$id,$mid,$seq,$con,$ac) = split(/\t/,$_);

   if($sp==9606) {
      $seed =~ s/U/T/g;
      $Aseq{$fam} = $seed;
      $Aprobe{$seq} = $fam;
   }
}
close(AFFY);
open (FAM,$miRcoords) or die;
while(<FAM>) {
   chomp();
   my($chr,$start,$end,$name,$score,$strand) = split(/\t/,$_);
   my($miR,$sequence) = split(":", $name);
   $seed = substr($sequence,1,7);
   $LU{$seed}{$miR} = 1;
}
close(FAM);

foreach $p(keys(%A0)) {
   if (exists($LU{$Aseq{$p}})) {
      @keys = keys (%{$LU{$Aseq{$p}}});
      foreach $k (@keys) {
	 print "$p\t$k\n";
      }

   }
   else {
      print "$p\tNF\n";
   }
}
