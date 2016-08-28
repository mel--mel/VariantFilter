#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <vector>
#include <map>
#include <type_traits>
#include <stdexcept> 

class  Variant{
public:
	static const int numElements = 8;

	Variant(std::string chrom, int pos, std::string id, std::string ref, std::string alt, 
			 float qual, std::string filter, std::map<std::string, std::string> info );
	~ Variant();

	static Variant string_to_variant(const std::string vcfLine);
	std::string variant_to_string();
	bool has_attribute(const std::string AttributeName, const std::string value, const std::string condition);
	bool has_attribute(const std::string AttributeName, const float Value, const std::string condition);


private:
	//assume 8 standard fields, of which the 8th field (info) has random number of attributes 
	std::string chrom;
	int pos;
	std::string id;
	std::string ref;
	std::string alt;
	float qual;
	std::string filter;
	//assume that each attribute is unique for each variant, so i can use a map to store them
	std::map<std::string, std::string> info; 
};

 Variant::Variant(std::string chrom, int pos, std::string id, std::string ref, std::string alt, 
			 float qual, std::string filter, std::map<std::string, std::string> info)
{
	this->chrom = chrom;
	this->pos = pos;
	this->id = id;
	this->ref = ref;
	this->alt = alt;
	this->qual = qual;
	this->filter = filter;
	this->info = info;
}

 Variant::~ Variant()
{
}

 Variant Variant::string_to_variant(const std::string vcfLine){
	std::string chrom;
	int pos;
	std::string id;
	std::string ref;
	std::string alt;
	float qual;
	std::string filter;
	std::map<std::string, std::string> info; 

	std::string helpstring;

	//get the first 8 fields of the line, separated by tabs
	//assume that the order of the fields is standard
	std::stringstream linestream(vcfLine);
	std::getline(linestream, chrom, '\t');
	//if an element is not a string, use a helping string and then convert it to the type needed
	std::getline(linestream, helpstring, '\t');
	pos = std::stoi(helpstring);
	std::getline(linestream, id, '\t');
	std::getline(linestream, ref, '\t');
	std::getline(linestream, alt, '\t');
	std::getline(linestream, helpstring, '\t');
	qual = std::stof(helpstring);
	std::getline(linestream, filter, '\t');
	
	
	std::getline(linestream, helpstring, '\t');
	//I now have the whole "info" as a string and i need to divide it into pairs (separated by ;)
	std::stringstream infostream(helpstring);
	std::string pairstring;
	std::size_t found;
	std::pair<std::string, std::string> pair;
	//ultil we reach the end of the infostream
	while(!((std::getline(infostream, pairstring, ';')).fail())){
		//divide each pair of string to pair.first and second - according to the = separator
		found = pairstring.find('=');
		if (found !=std::string::npos){
			pair.first = pairstring.substr(0, found);
			pair.second = pairstring.substr(++found, pairstring.length());
		}
		//if the given element is a flag, so it has not a value, fill pair.second with '.'
		else{
			pair.first = pairstring;
			pair.second = '.';
		}
		info.insert(pair);
	}

	//pass the elements required to the variant constructor
	return Variant(chrom, pos, id, ref, alt, qual, filter, info);
 }

 std::string Variant::variant_to_string(){
	 std::string variantstring;

	 //append the first seven "simple" fields to the string
	 variantstring.append(this->chrom);
	 variantstring.append("\t");
	 variantstring.append(std::to_string(this->pos));
	 variantstring.append("\t");
	 variantstring.append(this->id);
	 variantstring.append("\t");
	 variantstring.append(this->ref);
	 variantstring.append("\t");
	 variantstring.append(this->alt);
	 variantstring.append("\t");
	 variantstring.append(std::to_string(this->qual));
	 variantstring.append("\t");
	 variantstring.append(this->filter);
	 variantstring.append("\t");
	 
	 //append the "info" field of the variant, which should be constructed from the map
	 for (std::map<std::string, std::string>::iterator pair = this->info.begin(); pair != this->info.end(); ++pair){
		 //if in the middle of info, separate pairs by ";"
		 if (pair != this->info.begin()){
			variantstring.append(";");
		 }
		 //if this element of "info" is a flag, use nly pair.first (no value)
		 if (pair->second == "."){
			 variantstring.append(pair->first);
		 }
		 //else, it is an attribute name and value pair, separated by "="
		 else{
			 variantstring.append(pair->first);
			 variantstring.append("=");
			 variantstring.append(pair->second);
		 }
	 }

	 
	 return variantstring;
 }

 //compare floats
 bool is(const float actualValue, const float value, const std::string condition){
	if (condition.compare(">") == 0){
		 return (actualValue > value) ? true:false;
	}
	else if (condition.compare(">=") == 0){
		 return (actualValue >= value) ? true:false;
	}
	else if (condition.compare("<") == 0){
		return (actualValue < value) ? true:false;
	}
	else if (condition.compare("<=") == 0){
		return (actualValue <= value) ? true:false;
	}
	else if (condition.compare("==") == 0){
		return (actualValue == value) ? true:false;
	}
	else if (condition.compare("!=") == 0){
		return (actualValue != value) ? true:false;
	}
	else{
		throw "Invalid comparison condition! Conditions for floats should be: < or <= or > or >= or == or != .\r";
	}
 }

 //compare strings
 bool is(const std::string actualValue, const std::string value, const std::string condition){
	if (condition.compare(">") == 0){
		 return (actualValue.compare(value) > 0) ? true:false;
	}
	else if (condition.compare("<") == 0){
		return (actualValue.compare(value) < 0) ? true:false;
	}
	else if (condition.compare("==") == 0){
		return (actualValue.compare(value) == 0) ? true:false;
	}
	else if (condition.compare("!=") == 0){
		return (actualValue.compare(value) != 0) ? true:false;
	}
	else{
		throw "Invalid comparison condition! Conditions for strings should be: < or > or == or != .\r";
	}
 }

  bool Variant::has_attribute(const std::string attributeName, const float value, const std::string condition){
	 //search info map by attribute name
	 std::map<std::string, std::string>::iterator pair = this->info.find(attributeName);
	 //if attribute found
	 if (pair != this->info.end()){
		 try{
			 //check if its value can be converted to float
			float floatAttributeValue = std::stof(pair->second);
			//if yes, compare to given value
			return is(floatAttributeValue, value, condition);}
		 catch(const std::invalid_argument&){
			 //else throw exception
			 throw "Some attribute has a type that cannot be converted to float so filter comparisons are not valid!\r"; 
		 }
	 }
	 else{
		 return false;
	 }
 }

 bool Variant::has_attribute(const std::string attributeName, const std::string value, const std::string condition){
	 std::map<std::string, std::string>::iterator pair = this->info.find(attributeName);
	 if (pair != this->info.end()){
		 return is(pair->second, value, condition);
	 }
	 else{
		 return false;
	 }
 }

 class VariantCollection
 {
 public:
	 VariantCollection(std::string metainfo, std::string header, std::vector<Variant> variants);
	 ~VariantCollection();
	 static VariantCollection read_file(std::string FileName);
	 void write_file(std::string FileName);
	 void filter(std::string AttributeName, std::string condition, std::string AttributeValue);
	 void delete_variants(std::string AttributeName, std::string condition, std::string AttributeValue);
	 void delete_variants(std::string AttributeName, std::string condition, float floatAttributeValue);
	 
	 static bool is_it_metainfo(const std::string line);
	 static bool is_it_header(const std::string line);
	 std::string make_header(const std::string line);
	 std::string variants_to_string();

 private:
	 std::string MetaInfo;
	 std::string Header;
	 std::vector<Variant> Variants;

 };

 VariantCollection::VariantCollection(std::string metainfo, std::string header, std::vector<Variant> variants)
 {
	 MetaInfo = metainfo; 
	 Header = header;
	 Variants = variants;
 }

 VariantCollection::~VariantCollection()
 {
 }

 //each meta-info line starts with "##"
 bool VariantCollection::is_it_metainfo(std::string line){
	 return line.compare(0, 2, "##") == 0 ;
 }

 //the field header line starts with "#" 
 bool VariantCollection::is_it_header(const std::string line){
	 return ( line.compare(0, 1, "#") == 0 && line.compare(1, 2, "#") != 0 );
 }

 std::string make_header_line(const std::string line, const int numElements){
	 //keep the first "attributeNum" (= 8) fields of the header
	//the others are ignored, as stated in project requirements
	std::string header;
	 std::string line_to_search = line;
	std::size_t found = 0;
	//search for the first "attributeNum" elements of the header
	//if header has less elements, just copy the whole header
	for (int i = 0; i < numElements;){
		found = line_to_search.find('\t');
		if (found != std::string::npos){
			header.append(line_to_search.substr(0, ++found));
			line_to_search = line_to_search.substr(found, line_to_search.length());
			i++;
			}
		else{
			//if we reach the end of line earlier than expected, 
			//we copy the whole line to the header and terminate the loop
			header = line;
			i = numElements;
		}
					
	}
	return header;
 }

 VariantCollection  VariantCollection::read_file(const std::string FileName){
	 std::ifstream VcfFile;
	 std::string line;
	 std::string meta_info;
	 std::string header;
	 std::vector<Variant> variants;

	 VcfFile.open(FileName);
	if (VcfFile.is_open()){
		while ( std::getline (VcfFile,line) ){
			if (is_it_metainfo(line)){
				//I did not encode meta-information lines, time was not enough and I did not think
				//it was necessary for the project
				meta_info.append(line);
				meta_info.append("\r");
			}
			else if (is_it_header(line))
			{
				header.append(make_header_line(line, Variant::numElements));
				header.append("\r");
			}
			else{
				//if line is variant line, encode a variant and push it back into 
				//the variant vector of the class
				variants.push_back(Variant::string_to_variant(line));
			}
		}
		VcfFile.close();
		return VariantCollection(meta_info, header, variants);
	}
	else
	{
		throw "Could not open file for input!";
	}

 }

 std::string VariantCollection::variants_to_string(){
	 std::string variantstring;

	//for each variant in our variant collection
	for (std::vector<Variant>::iterator variant = this->Variants.begin(); variant != this->Variants.end(); ++variant){
		//get a line of this variant and append it to all variants string
		variantstring.append(variant->variant_to_string());
		variantstring.append("\r");
	}
	return variantstring;
 }

 void VariantCollection::write_file(std::string FileName){
	 std::ofstream NewVcfFile;

	 NewVcfFile.open(FileName);
	if (NewVcfFile.is_open()){
		//copy meta info and header strings as they are into output file
		NewVcfFile << this->MetaInfo;
		NewVcfFile << this->Header;
		//encode variant vector of the class into a single string and pass it to the output file
		NewVcfFile << variants_to_string();
		NewVcfFile.close();}
	else{
		throw "Could not open file for output!";}


 }

 void VariantCollection::delete_variants(std::string AttributeName, std::string condition, float floatAttributeValue){
	 //iterate through variants
	 for(std::vector<Variant>::iterator variant = this->Variants.begin(); variant != this->Variants.end();){
		 //if this variant has an attribute named like that, complying to this condition
		 if (variant->has_attribute(AttributeName, floatAttributeValue, condition)){
			 //delete this variant
			 variant = (this->Variants).erase(variant);
		 }
		 else //move to the next
		 {
			 ++variant;
		 }
	 }

 }

 void VariantCollection::delete_variants(std::string AttributeName, std::string condition, std::string AttributeValue){
	 //iterate through variants
	 for(std::vector<Variant>::iterator variant = this->Variants.begin(); variant != this->Variants.end();){
		 //if this variant has an attribute named like that, complying to this condition
		 if (variant->has_attribute(AttributeName, AttributeValue, condition)){
			 //delete this variant
			 variant = (this->Variants).erase(variant);
		 }
		 else //move to the next
		 {
			 ++variant;
		 }
	 }

 }

 void VariantCollection::filter(std::string AttributeName, std::string condition, std::string AttributeValue){
	 bool compare_to_float;
	 float floatAttributeValue;

	 //check if value to compare to attribute values is float
	 try{
		 floatAttributeValue = std::stof(AttributeValue);
		 compare_to_float = true;}
	 catch (const std::invalid_argument& ){
		 compare_to_float = false;}

	 //whether or not it is a float, call appropriate overload
	 if (compare_to_float){
		 this->delete_variants(AttributeName, condition, floatAttributeValue);
	 }
	 else{
		 this->delete_variants(AttributeName, condition, AttributeValue);
	 }

	 
 }

int main(){


	try{
		//assume that input file is encoded correctly
		VariantCollection VcfFile = VariantCollection::read_file("input.vcf");
		//filter out all variants with attribute AF > 0.5
		//by changing the attributes of this method we can change the filter
		//we can filter for every attribute we might want to, comparing to string or float values
		//with specific operands (> < >= <= == != for floats, > < == != for strings)
		//float operations are case sensitive
		VcfFile.filter("AF", ">", "0.5");
		//write filtered output to file
		VcfFile.write_file("newVcfFile.vcf");
	}
	catch (const char* exc){
		std::cout << exc << '\r';
	}

	return 0;
}