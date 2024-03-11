import sqlite3
import os

# Подключение к базе данных
conn = sqlite3.connect('patents.db')
cursor = conn.cursor()


# Функция для удаления записей по условию
def delete_invalid_compounds():
    # Получаем записи, где путь к файлу не "not"
    cursor.execute('SELECT id, structural_formula, molfile_path FROM compounds WHERE structural_formula != "not"')
    rows = cursor.fetchall()

    # Счетчик удаленных записей
    deleted_count = 0

    # Перебираем записи
    for row in rows:
        compound_id, structural_formula, molfile_path = row
        # Проверяем, существует ли файл по указанному пути
        if not os.path.exists(structural_formula):
            # Удаляем запись из таблицы compounds
            cursor.execute('DELETE FROM compounds WHERE id = ?', (compound_id,))
            # Удаляем связанные записи из таблицы compoundsInPatent
            cursor.execute('DELETE FROM compoundsInPatent WHERE compound_id = ?', (compound_id,))
            deleted_count += 1

    # Фиксируем изменения в базе данных
    conn.commit()

    return deleted_count


# Вызываем функцию и выводим количество удаленных записей
deleted_count = delete_invalid_compounds()
print(f'Удалено {deleted_count} соединений')

# Закрываем соединение с базой данных
conn.close()
