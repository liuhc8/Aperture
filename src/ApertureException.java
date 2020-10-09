
class ApertureException extends Exception{
	public ApertureException(String str) {
		super(str);
	}
}

class IllegalPositionException extends ApertureException{
	public IllegalPositionException(String str) {
		super(str);
	}
}

class IllegalBioFileException extends ApertureException{
	public IllegalBioFileException(String str) {
		super(str);
	}
}

class IllegalThreadNumberException extends ApertureException{
	public IllegalThreadNumberException(String str) {
		super(str);
	}
}

class IllegalBarcodeException extends ApertureException{
	public IllegalBarcodeException(String str) {
		super(str);
	}
}