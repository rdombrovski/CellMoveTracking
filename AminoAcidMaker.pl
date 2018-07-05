#!/usr/bin/perl

use strict; 
use warnings;

#amino acid molecular weights array
#Leucine and Isoleucine both have the same Molecular Weight so they have been combined into one amino acid result
my @MW = (89.0929, 133.1024, 147.1289, 165.1887, 75.0664, 155.1542, 131.1724, 
			146.1870, 132.1176, 115.1301, 146.1441, 174.2004, 105.0923, 119.1188, 
			117.1459, 204.2247, 181.1881);
		
#amino acid names array	
#Leucine and Isoleucine both have the same Molecular Weight so they have been combined into one amino acid result
my @AA = ("Alanine", "Aspartic Acid", "Glutamic Acid", "Phenylalanine", "Glycine", "Histidine", "Isoleucine/Leucine", "Lysine", "Aspargine",
			"Proline", "Glutamine", "Arginine", "Serine", "Threonine", "Valine", "Tryptophan", "Tyrosine");			

#arrays to store values for the closest 10 compounds to amino acids
my @closeList = ();
my @closeNums = ();
my @lastTenNums = ();
my @lastTenMol = ();

my @listMW = (1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000);
my @listComps = ();

#closest match is array including distance of compound weight from amino acid weight, compound name and amino acid name
my @closestMatch = (100, "c", "a");

#arrays for the 4 atom informations
my @atomMW = (1.0079, 15.9994, 14.0067, 12.0107);
my @currentAtom = ("H", "O", "N", "C");

my $isAA = 0;
my $numComps=0;


#loop until Amino Acid is formed
while ($isAA == 0){
	
my $numC=0;
my $numH=0;
my $numN=0;
my $numO=0;
my $compMW=0;
my $numAtom =0;
my $isComp = 0;
my $totalBV=0;
my $currentBV=0;
my $i=0;
my $compName="";
my $currentClose="";
my $closeW = 0;
my $j = 0;
my $smallDiff = 100;

print "so far, $numComps compounds have been made\n";

	#loop until compound is formed
	while ($isComp ==0 && $compMW < 204.2247){
		
		#pick a random number and action depending on which atom is chosen
		my $randNum = int(rand(4));
		
		if ($randNum==0){
			$numH+=1;
			$currentBV=1;
			$compMW = $compMW + $atomMW[0];
		}
		elsif ($randNum==1){
			$numO+=1;
			$currentBV=-2;
			$compMW = $compMW + $atomMW[1];
		}
		elsif ($randNum==2){
			$numN+=1;
			$currentBV=-3;
			$compMW = $compMW + $atomMW[2];
		}
		elsif ($randNum==3){
			$numC+=1;
			$currentBV=-4;
			$compMW = $compMW + $atomMW[3];
		}		
		
		$numAtom++;
		
		#add two to the Bond Value for each covalent bond
		if ($numAtom>1 && $currentBV<0 && $totalBV<0){
			$totalBV+=2;
		}	
				
		$totalBV += $currentBV;	
		
		#double bonds for every third Carbon
		if($numC != 0 && $numC % 3==0 && $randNum ==3){
			$totalBV+=2;
		}
			
		#when Bond Value gets below -4, add 3 hydrogens to the compound
		if ($totalBV < -4){
			$totalBV+=3;
			$numH +=3;
			$numAtom+=3;
			$compMW += 3 * $atomMW[0];			
		}
		
		#if bond value gets to zero, a compound is formed
		if($totalBV == 0){
			$isComp =1;
		}
		
		#if the compound's weight is above that of Tryptophan, then compound is formed
		if($compMW > 204.2247){
			$isComp =1;
		}

		#check if the compound is an amino acid, even when bond value is not 0 and 
		#compound does not have a MW above Tryptophan
		for ($i=0; $i<scalar @MW; $i++){
			if ($compMW == $MW[$i]){
				$isComp=1;
				$isAA =1;
				print "The first amino acid was : ", $AA[$i], " with a molar weight of ", $MW[$i], "\n";
			}	
		}
	}
	
	
	$numComps++;
	
	#acquire name of the compound
	if ($isComp==1){
		#print "The compound we had was: ";
		if ($numC!=0){
			$compName .= "C";
			$compName .= "$numC";
		}
		if ($numH!=0){
			$compName .= "H";
			$compName .= "$numH";		
		}
		if ($numN!=0){
			$compName .= "N";
			$compName .= "$numN";
		}
		if ($numO!=0){
			$compName .= "O";
			$compName .= "$numO";
		}
		#print $compName, "\n";
		#print "MW of compound was: ", $compMW, "\n";
	}
	
	#loop to acquire closest amino acid match to the compound
	for ($i=0;$i<scalar (@MW);$i++){
		my $diff = abs($compMW - $MW[$i]);
		if ($diff<=$smallDiff){
			$smallDiff = $diff;
			$currentClose = $AA[$i];
			$closeW = $MW[$i];
		}
	}
	
	#closest match so far for compound to amino acid
	if ($smallDiff<$closestMatch[0]){
		$closestMatch[0] = $smallDiff;
		$closestMatch[1] = $compName;
		$closestMatch[2] = $currentClose;
	}
	
	print "So far, the closest we got to an amino acid was $closestMatch[1], closing on $closestMatch[2] with a distance of $closestMatch[0]\n";
	
	
	
	
	#add a compound to the last ten list after each compound is formed
	push @lastTenMol, $compName;
	push @lastTenNums, $compMW;
	
	#every ten compounds formed, print out closest matches to amino acids to date
	if ($numComps %10 == 0){
			
			print "Every ten compounds summary!\n";
			#for the last ten compounds comparison loops
			for ($i=0; $i < scalar(@MW);$i++){
				my $smallTenDiff = 100;
				for ($j =0; $j < scalar(@lastTenNums); $j++){
					my $diff = abs($MW[$i] - $lastTenNums[$j]);
					
					if ($diff < $smallTenDiff){
						$smallTenDiff = $diff;
						$closeNums[$i] = $lastTenNums[$j];
						$closeList[$i] = $lastTenMol[$j];
					}
				}
			}
		
		#loop for comparing last ten compounds to all closest matches to date
		for ($i = 0; $i<scalar(@closeNums);$i++){
			
			my $diffNew = abs($MW[$i] - $closeNums[$i]);
			my $diffOld = abs($MW[$i] - $listMW[$i]);
			
			if($diffNew < $diffOld){
				$listMW [$i] = $closeNums[$i];
				$listComps[$i] = $closeList[$i]; 
			}
			
		}
	
	#print closest matches to all amino acids
		for ($i=0;$i<scalar(@closeNums);$i++){
			print "The closest compound to Amino Acid: $AA[$i] (MW=$MW[$i]) was $listComps[$i] (MW=$listMW[$i])\n";
		}
		
		
	#to check if amino acid is formed yet
	for ($i=0; $i<scalar (@MW); $i++){
		if ($compMW == $MW[$i]){
			$isAA =1;
			print "The first amino acid was : ", $AA[$i], " with a molar weight of ", $MW[$i], "\n";
		}	
	}
		
		@lastTenMol = [];
		@lastTenNums = [];	
	}

	
}


print "this is how many compounds it took : ", $numComps, "\n";
print "----------------------------------------------------------\n";