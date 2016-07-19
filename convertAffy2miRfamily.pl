#!/usr/bin/perl -w

$Affylist=shift;
open (LIST,$Affylist) or die;
while(<LIST>) {
   chomp();
   $A0{$_} =1;
}
close(LIST);

$AffyFile='Affy-mmu-miR-probes.txt';
$miRFamFile='miR_Family_Info.txt';
open (AFFY,$AffyFile) or die;
$head=<AFFY>;
while(<AFFY>) {
   chomp();
   my($probe,$id,$align,$len,$seq) = split(/\t/,$_);

   $Aseq{$probe} = substr($seq,1,7);
   $Aprobe{$seq} = $probe;
}
close(AFFY);
open (FAM,$miRFamFile) or die;
$head=<FAM>;
while(<FAM>) {
   chomp();

   my($fam,$seed,$sp,$id,$mid,$seq,$con,$ac) = split(/\t/,$_);

   if($sp==10090) {
      $LU{$seed} = $fam;
   }
}
close(FAM);

foreach $p(keys(%A0)) {
#print "$p\t";
   if (exists($LU{$Aseq{$p}})) {
      print "$LU{$Aseq{$p}}\n";
   }
   else {
      print "NF\n";
   }
}
