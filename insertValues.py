import sqlite3

# Подключаемся к базе данных
conn = sqlite3.connect('patents.db')
cursor = conn.cursor()

# Генерируем список molfile_path
molfile_paths = [f"/molfiles/US06175007/US06175008-20010116-C{i:05d}.MOL" for i in range(1, 228)]

# Формируем SQL запрос с использованием параметризации
sql = "INSERT INTO compounds (molfile_path) VALUES (?);"

# Вставляем строки в таблицу
cursor.executemany(sql, [(path,) for path in molfile_paths])

# Подтверждаем изменения и закрываем соединение
conn.commit()
conn.close()